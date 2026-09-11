#!/usr/bin/env python3
"""Compare the experimental Metal field updates with the SSE engine.

Run this script on macOS after installing an openEMS build configured with
-DWITH_METAL=ON.  It reports complete simulation time; this deliberately
includes command submission and synchronization, not just shader time.
"""

import argparse
import importlib.util
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time

# A source checkout keeps its Python dependencies in install/venv. Make the
# documented `python .../metal_fields.py` invocation use that interpreter when
# the caller's Python cannot import CSXCAD.
if importlib.util.find_spec('CSXCAD') is None:
    project_root = Path(__file__).resolve().parents[3]
    project_python = project_root / 'install' / 'venv' / 'bin' / 'python'
    if project_python.exists() and Path(sys.executable).resolve() != project_python.resolve():
        os.execv(str(project_python), [str(project_python), *sys.argv])
    raise ModuleNotFoundError(
        'CSXCAD is unavailable. Run with install/venv/bin/python or install '
        'the openEMS Python bindings for this interpreter.')

import h5py
import numpy as np

from CSXCAD import ContinuousStructure
from openEMS import openEMS


def make_model(path, cells, timesteps, nonuniform=False, boundaries=None, frequency=1e9):
    csx = ContinuousStructure()
    fdtd = openEMS(NrTS=timesteps, EndCriteria=0)
    fdtd.SetCSX(csx)
    fdtd.SetBoundaryCond(boundaries if boundaries is not None else ['PEC'] * 6)
    fdtd.SetGaussExcite(0, frequency)

    grid = csx.GetGrid()
    grid.SetDeltaUnit(1e-3)
    for axis, count in zip('xyz', cells):
        lines = np.arange(count + 1, dtype=float)
        if nonuniform:
            lines = count * (lines / count) ** 1.3
        grid.SetLines(axis, lines)

    excite = csx.AddExcitation('excite', 0, [1, 1, 1])
    center = np.asarray(cells, dtype=float) / 2
    excite.AddBox(center, center + 1)

    # Exercise nonuniform update coefficients, not only vacuum cells.
    dielectric = csx.AddMaterial('dielectric', epsilon=3.66, kappa=0.02)
    dielectric.AddBox(np.asarray(cells) * 0.2, np.asarray(cells) * 0.45)

    e_dump = csx.AddDump('Et', dump_type=0, dump_mode=0, file_type=1)
    e_dump.AddBox([0, 0, 0], list(cells))
    h_dump = csx.AddDump('Ht', dump_type=1, dump_mode=0, file_type=1)
    h_dump.AddBox([0, 0, 0], list(cells))
    fdtd.Write2XML(str(path))


def run(binary, model, engine, output, fp64_reference=False, compress=None, pml=None):
    output.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    if fp64_reference:
        env['OPENEMS_METAL_FP64_REFERENCE'] = '1'
    if compress is not None:
        env['OPENEMS_METAL_COMPRESS'] = '1' if compress else '0'
    if pml is not None:
        env['OPENEMS_METAL_PML'] = '1' if pml else '0'
    start = time.perf_counter()
    proc = subprocess.run(
        [binary, str(model), '--engine=' + engine], cwd=output, env=env,
        text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
    )
    elapsed = time.perf_counter() - start
    if proc.returncode:
        print(proc.stdout)
        raise RuntimeError('{} failed with exit code {}'.format(engine, proc.returncode))
    return elapsed, proc.stdout


def h5_arrays(path):
    values = {}
    with h5py.File(path, 'r') as h5:
        def collect(name, obj):
            if isinstance(obj, h5py.Dataset) and np.issubdtype(obj.dtype, np.number):
                values[name] = obj[...]
        h5.visititems(collect)
    return values


def compare(reference, result, rtol, atol):
    ref = h5_arrays(reference)
    got = h5_arrays(result)
    if ref.keys() != got.keys():
        raise AssertionError('HDF5 datasets differ: {} != {}'.format(ref.keys(), got.keys()))

    stats = {'max_abs': 0.0, 'max_rel': 0.0, 'relative_l2': 0.0,
             'rms': 0.0, 'samples': 0, 'failures': 0}
    squared_error = 0.0
    squared_reference = 0.0
    for name in ref:
        if ref[name].shape != got[name].shape:
            raise AssertionError('{} shape differs'.format(name))
        if not np.issubdtype(ref[name].dtype, np.floating):
            np.testing.assert_array_equal(got[name], ref[name], err_msg='dataset {}'.format(name))
            continue
        delta = np.abs(ref[name] - got[name])
        magnitude = np.abs(ref[name])
        tolerance = atol + rtol * magnitude
        # Pointwise relative error is meaningful only away from numerical zero.
        significant = magnitude > max(atol, float(magnitude.max(initial=0)) * 1e-6)
        stats['max_abs'] = max(stats['max_abs'], float(delta.max(initial=0)))
        if np.any(significant):
            stats['max_rel'] = max(stats['max_rel'], float(
                (delta[significant] / magnitude[significant]).max(initial=0)))
        squared_error += float(np.sum(delta.astype(np.float64) ** 2))
        squared_reference += float(np.sum(ref[name].astype(np.float64) ** 2))
        stats['samples'] += delta.size
        stats['failures'] += int(np.count_nonzero(delta > tolerance))
    stats['rms'] = np.sqrt(squared_error / stats['samples']) if stats['samples'] else 0.0
    stats['relative_l2'] = np.sqrt(squared_error / squared_reference) if squared_reference else 0.0
    return stats


def run_case(args, cells, timesteps, label):
    root = Path(tempfile.mkdtemp(prefix='openems-metal-fields-'))
    try:
        model = root / 'model.xml'
        make_model(model, cells, timesteps, args.nonuniform)
        sse_time, _ = run(args.openems, model, 'sse', root / 'sse')
        metal_time, metal_log = run(args.openems, model, 'metal', root / 'metal',
                                    args.fp64_reference,
                                    True if args.compare_dense else None)
        if args.compare_dense:
            if not any(message in metal_log for message in (
                    'Metal: lossless coefficients:',
                    'Metal: coefficient dictionary limit reached; using dense coefficients')):
                raise AssertionError('Coefficient compression/fallback was not exercised')
            _, dense_log = run(args.openems, model, 'metal', root / 'dense',
                               args.fp64_reference, compress=False)
            if 'Metal: lossless coefficients:' in dense_log:
                raise AssertionError('Dense Metal run unexpectedly enabled compression')
            for field in ('Et.h5', 'Ht.h5'):
                dense = h5_arrays(root / 'dense' / field)
                packed = h5_arrays(root / 'metal' / field)
                if dense.keys() != packed.keys():
                    raise AssertionError('Dense/compressed datasets differ')
                for name in dense:
                    a, b = dense[name], packed[name]
                    if a.shape != b.shape or a.dtype != b.dtype or a.tobytes() != b.tobytes():
                        raise AssertionError('Dense/compressed bits differ: ' + field + '/' + name)
            if args.fp64_reference:
                prefix = 'Metal FP64 update reference:'
                if ([s for s in dense_log.splitlines() if s.startswith(prefix)] !=
                        [s for s in metal_log.splitlines() if s.startswith(prefix)]):
                    raise AssertionError('Dense/compressed FP64 diagnostics differ')
        e = compare(root / 'sse' / 'Et.h5', root / 'metal' / 'Et.h5',
                    args.rtol, args.atol)
        h = compare(root / 'sse' / 'Ht.h5', root / 'metal' / 'Ht.h5',
                    args.rtol, args.atol)

        print('{}: {} x {} x {}, {} timesteps'.format(label, *cells, timesteps))
        if args.compare_dense:
            print('  Dense/compressed Metal: bit-identical complete E/H dumps')
            for line in metal_log.splitlines():
                if line.startswith(('Metal: lossless coefficients:',
                                    'Metal: coefficient dictionary limit')):
                    print('  ' + line)
        print('  SSE/Metal: {:.3f} / {:.3f} s ({:.2f}x)'.format(
              sse_time, metal_time, sse_time / metal_time))
        print('  E max abs/rel, relative L2, RMS, failures: '
              '{:.6g} / {:.6g}, {:.6g}, {:.6g}, {}/{}'.format(
              e['max_abs'], e['max_rel'], e['relative_l2'], e['rms'],
              e['failures'], e['samples']))
        print('  H max abs/rel, relative L2, RMS, failures: '
              '{:.6g} / {:.6g}, {:.6g}, {:.6g}, {}/{}'.format(
              h['max_abs'], h['max_rel'], h['relative_l2'], h['rms'],
              h['failures'], h['samples']))
        if 'Metal field updates' not in metal_log:
            raise AssertionError('Metal engine was not selected')
        if args.fp64_reference:
            reference_lines = [line for line in metal_log.splitlines()
                               if line.startswith('Metal FP64 update reference:')]
            if not reference_lines:
                raise AssertionError('FP64 reference result was not reported')
            print('  ' + reference_lines[-1])
        return e, h
    finally:
        if args.keep:
            print('  kept test data in {}'.format(root))
        else:
            shutil.rmtree(root)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--openems', default=os.environ.get('OPENEMS_BIN', 'openEMS'))
    parser.add_argument('--cells', type=int, nargs=3, default=(128, 128, 122))
    parser.add_argument('--timesteps', type=int, default=100)
    parser.add_argument('--suite', action='store_true', help='run several correctness cases')
    parser.add_argument('--fp64-reference', action='store_true',
                        help='run the diagnostic CPU FP64 update reference alongside Metal')
    parser.add_argument('--rtol', type=float, default=2e-4)
    parser.add_argument('--atol', type=float, default=1e-6)
    parser.add_argument('--keep', action='store_true')
    parser.add_argument('--compare-dense', action='store_true',
                        help='require bit-identical dumps from dense and compressed Metal')
    parser.add_argument('--nonuniform', action='store_true',
                        help='use varying mesh spacings to exercise dictionary fallback')
    args = parser.parse_args()

    if args.suite:
        cases = [
            ('tiny-boundaries', (8, 7, 6), 20),
            ('odd-dimensions', (17, 19, 21), 100),
            ('dielectric-medium', (48, 47, 46), 300),
            ('long-run', (24, 23, 22), 1000),
        ]
        results = [run_case(args, cells, steps, label)
                   for label, cells, steps in cases]
        print('suite maxima: E abs {:.6g}, H abs {:.6g}; tolerance failures: E {}, H {}'.format(
              max(result[0]['max_abs'] for result in results),
              max(result[1]['max_abs'] for result in results),
              sum(result[0]['failures'] for result in results),
              sum(result[1]['failures'] for result in results)))
    else:
        run_case(args, tuple(args.cells), args.timesteps, 'single')


if __name__ == '__main__':
    main()
