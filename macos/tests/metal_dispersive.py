#!/usr/bin/env python3
"""Verify the dispersive-material setup reuses the Metal geometry winners.

Mode '0' sets OPENEMS_METAL_PEC=0, forcing both the CPU PEC mapping and the
CSXCAD per-row fallback in Operator_Ext_LorentzMaterial. Mode 'default' uses the
Metal geometry pass (primal and dual winners) which the extension consumes.
Outputs must be bit-identical.
"""
import argparse
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import time

import numpy as np

from CSXCAD import ContinuousStructure
from CSXCAD.CSProperties import CSPropLorentzMaterial, CSPropDebyeMaterial
from openEMS import openEMS

from metal_fields import h5_arrays


def fixture(path, kind):
    cells = (20, 18, 16)
    csx = ContinuousStructure()
    fdtd = openEMS(NrTS=400, EndCriteria=0)
    fdtd.SetCSX(csx)
    fdtd.SetBoundaryCond(['PEC'] * 6)
    fdtd.SetGaussExcite(0, 1e9)
    grid = csx.GetGrid()
    grid.SetDeltaUnit(1e-3)
    for axis, count in zip('xyz', cells):
        grid.SetLines(axis, np.arange(count + 1, dtype=float))
    excite = csx.AddExcitation('excite', 0, [1, 1, 1])
    excite.AddBox([10, 9, 8], [11, 10, 9])

    pset = csx.GetParameterSet()
    if kind == 'lorentz':
        # Electric and magnetic poles exercise both the primal and dual winners.
        mat = CSPropLorentzMaterial(pset, order=1, epsilon=2.5)
        mat.SetDispersiveMaterialProperty(
            0, eps_plasma=1e8, eps_pole_freq=2e8, eps_relax=1e-9,
            mue_plasma=5e7, mue_pole_freq=1e8, mue_relax=2e-9)
    else:
        mat = CSPropDebyeMaterial(pset, order=1, epsilon=3.0)
        mat.SetDispersiveMaterialProperty(0, eps_delta=0.5, eps_relax=3e-10)
    mat.SetName(kind)
    mat.AddBox([2, 2, 2], [12, 12, 12], priority=3)
    csx.AddProperty(mat)

    e_dump = csx.AddDump('Et', dump_type=0, dump_mode=0, file_type=1)
    e_dump.AddBox([0, 0, 0], list(cells))
    h_dump = csx.AddDump('Ht', dump_type=1, dump_mode=0, file_type=1)
    h_dump.AddBox([0, 0, 0], list(cells))
    fdtd.Write2XML(str(path))


def run(binary, model, mode, output):
    output.mkdir()
    env = os.environ.copy()
    if mode == 'default':
        env.pop('OPENEMS_METAL_PEC', None)
    else:
        env['OPENEMS_METAL_PEC'] = mode
    start = time.perf_counter()
    p = subprocess.run([binary, str(model), '--engine=metal'], cwd=output,
                       env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    wall = time.perf_counter() - start
    (output / 'solver.log').write_text(p.stdout)
    if p.returncode:
        raise RuntimeError(p.stdout)
    fast = 'resolved geometry winners' in p.stdout
    if mode == 'default':
        if 'Metal PEC:' not in p.stdout:
            raise AssertionError('Metal geometry pass was not exercised')
        if not fast:
            raise AssertionError('dispersive material did not use resolved winners')
        primal, dual = (int(v) for v in re.search(
            r'Metal dispersive material: (\d+) primal, (\d+) dual cells', p.stdout).groups())
        if dual == 0:
            raise AssertionError('dual (magnetic) winners were not resolved')
    elif fast:
        raise AssertionError('fallback mode unexpectedly used resolved winners')
    timings = [s for s in p.stdout.splitlines()
               if s.startswith(('Metal PEC:', 'CPU PEC:', 'Time for ', 'Metal dispersive material:'))]
    print(mode, 'wall {:.3f}s'.format(wall), '; '.join(timings), flush=True)
    return p.stdout


def compare(a, b):
    files = {p.relative_to(a) for p in a.rglob('*') if p.is_file() and p.name != 'solver.log'}
    other = {p.relative_to(b) for p in b.rglob('*') if p.is_file() and p.name != 'solver.log'}
    if files != other:
        raise AssertionError('Output file sets differ')
    for name in files:
        if name.suffix == '.h5':
            x, y = h5_arrays(a / name), h5_arrays(b / name)
            if x.keys() != y.keys():
                raise AssertionError('HDF5 datasets differ')
            for key in x:
                if (x[key].dtype != y[key].dtype or x[key].shape != y[key].shape
                        or x[key].tobytes() != y[key].tobytes()):
                    raise AssertionError(str(name) + '/' + key)
        else:
            rows = lambda p: [s for s in p.read_text().splitlines()
                              if s.strip() and not s.startswith(('#', '%'))]
            if rows(a / name) != rows(b / name):
                raise AssertionError(str(name))
    print('  Bit-identical field datasets / probe numeric outputs')


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--openems', required=True)
    p.add_argument('--keep', action='store_true')
    args = p.parse_args()
    root = Path(tempfile.mkdtemp(prefix='openems-metal-disp-'))
    try:
        for kind in ('lorentz', 'debye'):
            model = root / (kind + '.xml')
            fixture(model, kind)
            print(kind, flush=True)
            for mode in ('0', 'default'):
                run(str(Path(args.openems).resolve()), model, mode, root / (kind + '-' + mode))
            compare(root / (kind + '-0'), root / (kind + '-default'))
    finally:
        if args.keep:
            print('Kept', root)
        else:
            shutil.rmtree(root)


if __name__ == '__main__':
    main()
