# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 openEMS contributors
"""Native integration control: repeated immutable snapshots versus stopped prefixes.

Each engine is compared against itself: final complex sums and probe values must
agree with/without requests; checkpoint sums must match a run stopped at
that timestep within the tolerances below. Engine-to-engine rounding is
tested separately by GPU_Engine.py.
Run with --binary PATH --engine basic|gpu-reference|gpu --output DIR.
"""
import argparse
import json
import subprocess
import time
import xml.etree.ElementTree as ET
from pathlib import Path
import h5py
import numpy as np
from CSXCAD import ContinuousStructure
from openEMS import openEMS
from openEMS.snapshots import request_dft_snapshot


# These are same-engine controls with no intended arithmetic changes. Allow
# about eight float32 epsilons in the DFT sums for platform/compiler rounding;
# scale the absolute floor to each reference array's peak, not to SI units.
FIELD_RTOL = 1e-6
FIELD_ATOL_PEAK = 1e-7


def model(path, steps):
    fdtd = openEMS(NrTS=steps, EndCriteria=0)
    fdtd.SetGaussExcite(3e9, 2e9)
    fdtd.SetBoundaryCond(['PEC'] * 6)
    csx = ContinuousStructure()
    fdtd.SetCSX(csx)
    mesh = csx.GetGrid()
    mesh.SetDeltaUnit(1e-3)
    for ax in 'xyz':
        mesh.AddLine(ax, np.linspace(-30, 30, 25))
    source = csx.AddExcitation('source', exc_type=0, exc_val=[0, 0, 1])
    source.AddBox([0, 0, -2.5], [0, 0, 2.5])
    csx.AddProbe('et', p_type=2).AddPoint([5, 5, 0])
    line = [-20, 20]
    dumps = [
        csx.AddDump('nf2ff_' + field, dump_type=dump, dump_mode=1,
                    file_type=1, frequency=[3e9, 5e9])
        for field, dump in [('E', 10), ('H', 11)]
    ]
    for n in range(6):
        a, b = [line[0]] * 3, [line[1]] * 3
        a[n // 2] = b[n // 2] = line[n % 2]
        for dump in dumps:
            dump.AddBox(a, b)
    fdtd.Write2XML(str(path))


def fields(path):
    with h5py.File(path, 'r') as f:
        return {name: value[()] for name, value in f['FieldData/FD'].items()}


def compare(a, b):
    worst = 0.0
    nonzero = 0
    worst_relative = 0.0
    for field in 'EH':
        for face in range(6):
            name = 'nf2ff_%s_%d.h5' % (field, face)
            x, y = fields(a / name), fields(b / name)
            assert x.keys() == y.keys()
            for key in x:
                assert np.isfinite(x[key]).all()
                peak = float(np.max(np.abs(y[key])))
                error = float(np.max(np.abs(x[key] - y[key])))
                np.testing.assert_allclose(
                    x[key], y[key], rtol=FIELD_RTOL,
                    atol=FIELD_ATOL_PEAK * peak,
                    err_msg='%s/%s: max error=%g, reference peak=%g' %
                            (name, key, error, peak))
                worst = max(worst, error)
                if peak > 0:
                    worst_relative = max(worst_relative, error / peak)
                nonzero += int(np.any(np.abs(x[key]) > 0))
    return dict(max_absolute_error=worst, max_peak_relative_error=worst_relative,
                nonzero_arrays=nonzero)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--binary', required=True)
    parser.add_argument('--engine', default='basic')
    parser.add_argument('--output', required=True)
    parser.add_argument('--steps', type=int, default=30000)
    parser.add_argument('--legacy', action='store_true')
    parser.add_argument(
        '--nf2ff', action='store_true',
        help='run the existing CPU NF2FF tool while FDTD continues')
    args = parser.parse_args()
    base = Path(args.output).resolve()
    base.mkdir(parents=True, exist_ok=False)
    command = [str(Path(args.binary).resolve()), 'model.xml', '--engine=' + args.engine,
               '--numThreads=2', '--exact-endcriteria']
    if args.legacy:
        command += ['--legacyHDF5Dumps']
    live, control = base / 'live', base / 'control'
    live.mkdir()
    control.mkdir()
    model(live / 'model.xml', args.steps)
    (control / 'model.xml').write_bytes((live / 'model.xml').read_bytes())
    (control / 'DFT_SNAPSHOT').write_text('ignored\n')
    request_dft_snapshot(live, 'first')
    started = time.monotonic()
    preview = None
    with (live / 'solver.log').open('w') as log:
        process = subprocess.Popen(
            command + ['--dft-snapshots'], cwd=live,
            stdout=log, stderr=subprocess.STDOUT)
        try:
            for label in ['first', 'second']:
                if label == 'second':
                    while (live / 'DFT_SNAPSHOT.active').exists():
                        assert process.poll() is None
                        time.sleep(0.001)
                    request_dft_snapshot(live, label)
                checkpoint = live / 'dft-snapshots' / label
                while not (checkpoint / 'snapshot.xml').exists():
                    log_tail = (live / 'solver.log').read_text()[-3000:]
                    assert process.poll() is None, log_tail
                    assert not checkpoint.with_name(label + '.error').exists(), log_tail
                    assert time.monotonic() - started < 120
                    time.sleep(0.01)
                assert len(list(checkpoint.glob('*.h5'))) == 12
                # The published directory must be readable immediately, without SWMR.
                fields(checkpoint / 'nf2ff_E_0.h5')
                if label == 'second' and args.nf2ff:
                    # The transform is a different process; the FDTD worker is
                    # still running and is not paused by this consumer.
                    assert process.poll() is None
                    nfxml = ET.Element(
                        'nf2ff', freq='3000000000,5000000000', Center='0,0,0',
                        Outfile=str(base / 'live-preview.h5'), NumThreads='2')
                    ET.SubElement(nfxml, 'theta').text = ','.join(
                        map(str, np.deg2rad(np.arange(0, 181, 15))))
                    ET.SubElement(nfxml, 'phi').text = ','.join(
                        map(str, np.deg2rad(np.arange(0, 360, 30))))
                    for face in range(6):
                        ET.SubElement(
                            nfxml, 'Planes',
                            E_Field=str(checkpoint / ('nf2ff_E_%d.h5' % face)),
                            H_Field=str(checkpoint / ('nf2ff_H_%d.h5' % face)))
                    ET.ElementTree(nfxml).write(base / 'live-preview.xml')
                    before = time.monotonic()
                    with (base / 'preview.log').open('w') as nf_log:
                        subprocess.run(
                            [str(Path(args.binary).resolve().with_name('nf2ff')),
                             str(base / 'live-preview.xml')],
                            stdout=nf_log, stderr=subprocess.STDOUT,
                            timeout=60, check=True)
                    with h5py.File(base / 'live-preview.h5', 'r') as h:
                        assert np.isfinite(h['nf2ff/E_theta/FD/f0'][()]).all()
                    preview = dict(
                        backend='cpu', seconds=time.monotonic() - before,
                        fdtd_alive_at_start=True,
                        fdtd_alive_at_end=process.poll() is None)
            assert process.wait(timeout=120) == 0
        finally:
            if process.poll() is None:
                process.kill()
                process.wait()
    with (control / 'solver.log').open('w') as log:
        subprocess.run(command, cwd=control, stdout=log,
                       stderr=subprocess.STDOUT, check=True, timeout=120)
    assert (control / 'DFT_SNAPSHOT').read_text() == 'ignored\n'
    assert not (control / 'dft-snapshots').exists()
    checks = {'final': compare(live, control)}
    assert checks['final']['nonzero_arrays'] > 0
    # Probe output is float64 text; use the same relative allowance and no
    # absolute floor so tiny time values are not hidden by a field-unit floor.
    np.testing.assert_allclose(
        np.loadtxt(live / 'et', comments='%'),
        np.loadtxt(control / 'et', comments='%'),
        rtol=FIELD_RTOL, atol=0, err_msg='final probe output')
    steps = []
    for label in ['first', 'second']:
        checkpoint = live / 'dft-snapshots' / label
        header = ET.parse(checkpoint / 'snapshot.xml').getroot()
        assert header.get('converged') == 'false'
        step = int(header.get('timestep'))
        steps.append(step)
        prefix = base / ('prefix-' + label)
        prefix.mkdir()
        xml = ET.parse(live / 'model.xml')
        xml.getroot().find('FDTD').set('NumberOfTimesteps', str(step))
        xml.write(prefix / 'model.xml')
        with (prefix / 'solver.log').open('w') as log:
            subprocess.run(command, cwd=prefix, stdout=log,
                           stderr=subprocess.STDOUT, check=True, timeout=120)
        checks[label] = compare(checkpoint, prefix)
        with h5py.File(checkpoint / 'nf2ff_E_0.h5', 'r') as f:
            assert int(f.attrs['last_sample_timestep']) <= step
            assert int(f.attrs['sample_count']) > 0
    assert 0 < steps[0] < steps[1] < args.steps
    # Failure and traversal requests must not stop stepping or publish a set.
    for label, request_id in [('no-fields', 'empty'), ('invalid-id', '../escape')]:
        folder = base / label
        folder.mkdir()
        xml = ET.parse(live / 'model.xml')
        xml.getroot().find('FDTD').set('NumberOfTimesteps', '100')
        props = xml.getroot().find('ContinuousStructure/Properties')
        for dump in list(props.findall('DumpBox')):
            props.remove(dump)
        xml.write(folder / 'model.xml')
        (folder / 'DFT_SNAPSHOT').write_text(request_id + '\n')
        with (folder / 'solver.log').open('w') as log:
            subprocess.run(
                command + ['--dft-snapshots'], cwd=folder, stdout=log,
                stderr=subprocess.STDOUT, check=True, timeout=60)
        assert (folder / 'et').exists()
        assert not (folder / 'DFT_SNAPSHOT').exists()
        if label == 'no-fields':
            assert (folder / 'dft-snapshots/empty.error').exists()
            assert not (folder / 'dft-snapshots/empty').exists()
        else:
            assert not (folder / 'dft-snapshots').exists()
    result = dict(
        engine=args.engine, legacy=args.legacy, snapshot_steps=steps, checks=checks,
        field_rtol=FIELD_RTOL, field_atol_peak=FIELD_ATOL_PEAK,
        elapsed_seconds=time.monotonic() - started, live_nf2ff=preview)
    (base / 'result.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
