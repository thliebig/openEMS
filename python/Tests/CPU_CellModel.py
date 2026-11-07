# -*- coding: utf-8 -*-
"""
 CPU engine: the cell model (CellModel) of three structures with known numbers

 The CPU engine (engine='cpu') turns the operator into runs of cell classes
 along its fast axis. A mistake there need not change any result (Engine_Compare.py
 compares the fields): a lossless cell taken as lossy, runs broken where nothing
 changes, inert cells no longer skipped, UPML layers not merged into a profile
 all compute the same fields, only slower. This test checks the model itself,
 by the counts the engine prints (verbose), on structures whose numbers follow
 by hand:

   vacuum in a PEC box, strongly graded mesh: 9 coefficient sets, ZERO and
     LOSSLESS per voltage component, LOSSLESS per current component (PEC
     zeroes the voltages only); m = vi/g does not depend on the geometry in a
     homogeneous material, so the grading must not split sets
   a PEC block of L cells per side in vacuum: 3*(L+1)*L^2 inert currents; a
     current is inert where the four voltages of its curl lie in the metal,
     its surface edges included
   vacuum with a PML on the two z faces only, z the fast axis: the UPML layers
     along a line merge into one profile run per side, about 3 runs per line;
     without the profiles a run per layer, a mean run length below 2

 Run with --engine=cpu -v; `--debug-operator` writes the whole model to
 cellmodel_dump.h5 (CellModel::Dump(): class, set and run per edge, the sets,
 runs and poles per component) for a closer look, e.g. with h5py.

 Tested with
  - python 3.13
  - openEMS v0.37+

 (c) 2026 Thorsten Liebig <thorsten.liebig@gmx.de>
"""

import os, re, sys, tempfile, ctypes, ctypes.util
import numpy as np

from CSXCAD import ContinuousStructure
from openEMS import openEMS

N = 31   # mesh lines per axis


def model(lines, boundary, metal_block=None):
    FDTD = openEMS(NrTS=10, EndCriteria=0)
    FDTD.SetGaussExcite(5e9, 4e9)
    FDTD.SetBoundaryCond(boundary)
    CSX = ContinuousStructure()
    FDTD.SetCSX(CSX)
    mesh = CSX.GetGrid()
    mesh.SetDeltaUnit(1e-3)
    for d in 'xyz':
        mesh.AddLine(d, lines)
    if metal_block is not None:
        CSX.AddMetal('block').AddBox(*metal_block)
    # a short source away from the block
    CSX.AddExcitation('e', exc_type=0, exc_val=[0, 0, 1]).AddBox([lines[3]]*2 + [lines[2]], [lines[3]]*2 + [lines[4]])
    return FDTD


def counts(FDTD, name, env=None):
    """ the counts line of the CPU engine for a run of FDTD """
    Sim_Path = os.path.join(tempfile.gettempdir(), 'CPU_CellModel_' + name)
    log_file = Sim_Path + '.log'
    for k, v in (env or {}).items():
        os.environ[k] = v
    sys.stdout.flush()
    saved = os.dup(1)
    try:
        with open(log_file, 'w') as log:
            os.dup2(log.fileno(), 1)
            try:
                FDTD.Run(Sim_Path, cleanup=True, engine='cpu', exact_endcriteria=True, verbose=1)
            finally:
                try:   # flush the C stdio buffer of the openEMS library
                    ctypes.CDLL(ctypes.util.find_library('c')).fflush(None)
                except (OSError, AttributeError, TypeError):
                    pass
                os.dup2(saved, 1)
    finally:
        os.close(saved)
        for k in (env or {}):
            del os.environ[k]
    with open(log_file) as log:
        text = log.read()
    m = re.search(r'^CellModel: counts: (.*)$', text, re.M)
    assert m, 'FAIL [{}]: no counts line of the cell model in the output of the cpu engine'.format(name)
    c = {k: int(v) for k, v in re.findall(r'(\w+) (\d+)', m.group(1))}
    print('{}: {}'.format(name, m.group(1)))
    return c


uniform = np.arange(N, dtype=float)

# 1) vacuum, PEC walls, graded mesh (cells from 0.25 to 4 mm): ZERO and LOSSLESS per voltage
#    component, LOSSLESS per current component
graded = np.cumsum(np.concatenate([[0], np.geomspace(0.25, 4, N-1)]))
c = counts(model(graded, ['PEC']*6), 'graded')
assert c['sets'] == 9, 'FAIL [graded]: {} coefficient sets, expected 9 (voltages ZERO and LOSSLESS, currents LOSSLESS)'.format(c['sets'])
assert c['lossy'] == 0 and c['upml'] == 0 and c['dispersive'] == 0, 'FAIL [graded]: classes other than ZERO and LOSSLESS: {}'.format(c)

# 2) a PEC block of L cells per side: 3*(L+1)*L^2 inert currents
L = 10
c = counts(model(uniform, ['PEC']*6, metal_block=([10, 10, 10], [10+L, 10+L, 10+L])), 'pec_block')
expected = 3*(L+1)*L*L
assert c['inert'] == expected, 'FAIL [pec_block]: {} inert currents, expected 3*(L+1)*L^2 = {}'.format(c['inert'], expected)

# 3) PML on the z faces, z the fast axis: the UPML layers merge into profiles
c = counts(model(uniform, ['PEC', 'PEC', 'PEC', 'PEC', 'PML_8', 'PML_8']), 'pml_profile', env={'OPENEMS_CPU_FAST_AXIS': 'z'})
cells = sum(c[k] for k in ('zero', 'lossless', 'lossy', 'upml', 'dispersive', 'inert'))
mean = cells / c['runs']
# about 3 runs per line of 31 cells (profile, vacuum, profile): ~10; without profiles a run
# per UPML layer, below 2. 6 leaves room for the lines at the walls.
assert mean > 6, 'FAIL [pml_profile]: mean run length {:.2f}, expected > 6 (the UPML layers merged into profiles)'.format(mean)

print('PASS')
