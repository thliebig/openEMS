# -*- coding: utf-8 -*-
"""
 Extension Optimizations Test

 The multithreaded engine skips the empty extension hooks and the SSE engines
 update the UPML through a cursor and row pointers. All of it is meant to be
 exact: the same setup is run with and without `--no-ext-opt` and the port
 voltage and current must match bit for bit. The setup has a UPML, a Mur
 boundary, a lumped port and a lumped RLC element, so every optimized extension
 type takes part.

 Pass criteria:
   port voltage and current files identical, with and without optimizations

 The comparison is exact on purpose: no tolerance is needed because the
 optimizations do not change the order of any floating point operation.
 If this test fails after a change to an extension or the engine, the
 optimization (or the change) altered the arithmetic.

 (c) 2026 Ismail Akdag

"""

import os, tempfile, shutil
import numpy as np

from CSXCAD  import ContinuousStructure
from openEMS import openEMS

Sim_Path = os.path.join(tempfile.gettempdir(), 'Extension_Optimizations')

def run(sub_path, **options):
    FDTD = openEMS(NrTS=3000, EndCriteria=0)
    FDTD.SetGaussExcite(2e9, 1e9)
    FDTD.SetBoundaryCond(['PML_8', 'PML_8', 'MUR', 'MUR', 'PML_8', 'PML_8'])

    CSX = ContinuousStructure()
    FDTD.SetCSX(CSX)
    mesh = CSX.GetGrid()
    mesh.SetDeltaUnit(1e-3)
    for d in 'xyz':
        mesh.AddLine(d, np.linspace(-12, 12, 25))

    # a lumped port with a wire dipole arm and a parallel RLC load
    port = FDTD.AddLumpedPort(1, 50, [0, 0, -1], [0, 0, 1], 'z', 1.0, priority=5)
    wire = CSX.AddMetal('wire')
    wire.AddBox([0, 0, 1], [0, 0, 6], priority=10)
    load = CSX.AddLumpedElement('load', ny='z', caps=False, R=100, L=20e-9, C=0.5e-12, LEtype=0)
    load.AddBox([-1, 0, 6], [1, 0, 7])

    path = Sim_Path + '_' + sub_path
    FDTD.Run(path, cleanup=True, numThreads=2, **options)
    return path

def port_data(path, name):
    # skip the comment header, which carries the date
    with open(os.path.join(path, name)) as f:
        return [line for line in f if not line.startswith('%')]

path_opt   = run('optimized')
path_plain = run('plain', no_ext_opt=True)

for name in ('port_ut_1', 'port_it_1'):
    opt   = port_data(path_opt, name)
    plain = port_data(path_plain, name)
    assert len(opt) > 50, 'FAIL: {} has only {} samples'.format(name, len(opt))
    assert opt == plain, 'FAIL: {} differs with --no-ext-opt'.format(name)

print('port voltage and current identical with and without the optimizations')
# Run() leaves the working directory inside the simulation folder
os.chdir(tempfile.gettempdir())
shutil.rmtree(path_opt)
shutil.rmtree(path_plain)
