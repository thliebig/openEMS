# -*- coding: utf-8 -*-
"""
 GPU engine test — bit-exact comparison with the basic engine

 The GPU engine (engine='gpu') currently runs the reference backend, which
 performs the main FDTD updates on the CPU on its own copy of the fields,
 exactly like the basic engine. Extensions without a device implementation run
 on the host copy of the fields, which the GPU engine synchronizes around every
 half-step. The results must therefore be bit-identical to the basic engine.

 Each case runs with engine='basic' and engine='gpu' and compares all probe
 files and all HDF5 field dumps. Together the cases use the excitation, UPML,
 Mur ABC, Lorentz material, lumped RLC, conducting sheet, TF/SF and
 steady-state extensions.

 Pass criteria (per case)
   the GPU engine was created for engine='gpu' (not a silent fallback)
   all probe data and field dumps are bit-identical to the basic engine

 Tested with
  - python 3.13
  - openEMS v0.0.37+

 (c) 2026 Sean Mollet <sean@malmoset.com>

"""

import os, sys, glob, tempfile, ctypes, ctypes.util
import numpy as np
import h5py

from CSXCAD  import ContinuousStructure
from CSXCAD.CSProperties import CSPropLorentzMaterial
from openEMS import openEMS
from openEMS.ports import LumpedPort

unit = 1e-3   # drawing unit: mm


def channel_1d(FDTD, CSX, sinus=False):
    """ TEM channel (PEC walls normal to x, PMC walls normal to y, PML at the z-ends) """
    if sinus:
        FDTD.SetSinusExcite(3e9)
    else:
        FDTD.SetGaussExcite(5.5e9, 4.5e9)
    FDTD.SetBoundaryCond(['PEC', 'PEC', 'PMC', 'PMC', 'PML_8', 'PML_8'])
    mesh = CSX.GetGrid()
    mesh.AddLine('x', [0, 0.5, 1])
    mesh.AddLine('y', [0, 0.5, 1])
    mesh.AddLine('z', np.arange(0, 100.5, 0.5))
    CSX.AddExcitation('plane', exc_type=0, exc_val=[1, 0, 0]).AddBox([0, 0, 20], [1, 1, 20])
    CSX.AddProbe('et', p_type=2).AddPoint([0.5, 0.5, 70])
    CSX.AddProbe('ht', p_type=3).AddPoint([0.5, 0.5, 70])


def case_dispersive_pml():
    FDTD = openEMS(NrTS=5000, EndCriteria=1e-5)
    CSX = ContinuousStructure()
    FDTD.SetCSX(CSX)
    CSX.GetGrid().SetDeltaUnit(unit)
    channel_1d(FDTD, CSX)
    drude = CSPropLorentzMaterial(CSX.GetParameterSet(), order=1)
    drude.SetName('drude')
    drude.SetDispersiveMaterialProperty(0, eps_plasma=5e9, eps_relax=1e-9)
    CSX.AddProperty(drude)
    drude.AddBox([0, 0, 40], [1, 1, 50], priority=10)
    return FDTD, CSX


def case_3d_mixed():
    FDTD = openEMS(NrTS=800, EndCriteria=0)
    FDTD.SetGaussExcite(5.5e9, 4.5e9)
    FDTD.SetBoundaryCond(['MUR'] * 6)
    CSX = ContinuousStructure()
    FDTD.SetCSX(CSX)
    mesh = CSX.GetGrid()
    mesh.SetDeltaUnit(unit)
    for ax in 'xyz':
        mesh.AddLine(ax, np.arange(-10, 10.5, 1))

    # lumped port and series RLC element connected by two wires
    LumpedPort(CSX, 1, 50, [-4, 0, 0], [-4, 0, 2], 'z', excite=1)
    rlc = CSX.AddLumpedElement('rlc', ny='z', caps=False, R=10, L=1e-9, C=1e-12, LEtype=1)
    rlc.AddBox([4, 0, 0], [4, 0, 2], priority=10)
    wire = CSX.AddMetal('wire')
    wire.AddCurve([[-4, 4], [0, 0], [0, 0]])
    wire.AddCurve([[-4, 4], [0, 0], [2, 2]])

    # conducting sheet patch
    sheet = CSX.AddConductingSheet('sheet', conductivity=5.8e7, thickness=1e-6)
    sheet.AddBox([-3, -3, 5], [3, 3, 5], priority=10)

    # TF/SF plane wave
    pw = CSX.AddExcitation('plane_wave', exc_type=10, exc_val=[0, 0, 1])
    pw.SetPropagationDir([1, 0, 0])
    pw.SetFrequency(5e9)
    pw.AddBox([-6, -6, -6], [6, 6, 6])

    CSX.AddProbe('et', p_type=2).AddPoint([0, 5, 0])
    CSX.AddProbe('ht', p_type=3).AddPoint([0, 5, 0])
    CSX.AddDump('Et', dump_type=0, file_type=1).AddBox([-10, -10, 0], [10, 10, 0])
    CSX.AddDump('Hf', dump_type=11, file_type=1, frequency=[5e9]).AddBox([-10, 0, -10], [10, 0, 10])
    return FDTD, CSX


def case_steady_state():
    FDTD = openEMS(NrTS=100000, EndCriteria=1e-6)
    CSX = ContinuousStructure()
    FDTD.SetCSX(CSX)
    CSX.GetGrid().SetDeltaUnit(unit)
    channel_1d(FDTD, CSX, sinus=True)
    return FDTD, CSX


def run_captured(case, Sim_Path, engine):
    """ Run a case with the given engine, return the openEMS console output """
    FDTD, CSX = case()
    log_file = Sim_Path + '.log'
    sys.stdout.flush()
    saved = os.dup(1)
    with open(log_file, 'w') as log:
        os.dup2(log.fileno(), 1)
        try:
            FDTD.Run(Sim_Path, cleanup=True, engine=engine)
        finally:
            try:   # flush the C stdio buffer of the openEMS library
                ctypes.CDLL(ctypes.util.find_library('c')).fflush(None)
            except (OSError, AttributeError, TypeError):
                pass
            os.dup2(saved, 1)
            os.close(saved)
    with open(log_file) as log:
        return log.read()


def compare_h5(fn_a, fn_b):
    """ return the names of all datasets that differ between two HDF5 files """
    diff = []
    with h5py.File(fn_a, 'r') as a, h5py.File(fn_b, 'r') as b:
        def visit(name, obj):
            if isinstance(obj, h5py.Dataset):
                if name not in b or not np.array_equal(obj[()], b[name][()]):
                    diff.append(name)
        a.visititems(visit)
    return diff


def compare_outputs(path_a, path_b):
    """ compare all probe files and HDF5 dumps, return a list of differences """
    diff = []
    files = [f for f in os.listdir(path_a) if os.path.isfile(os.path.join(path_a, f))]
    probes = [f for f in files if '.' not in f]   # probe and port files have no extension
    dumps  = [f for f in files if f.endswith('.h5')]
    assert probes, f'FAIL: no probe files in {path_a}'
    for f in probes:
        a = np.loadtxt(os.path.join(path_a, f), comments='%')
        b = np.loadtxt(os.path.join(path_b, f), comments='%')
        if not np.array_equal(a, b):
            diff.append(f)
    for f in dumps:
        diff += [f'{f}:{d}' for d in compare_h5(os.path.join(path_a, f), os.path.join(path_b, f))]
    return diff, len(probes), len(dumps)


cases = [('dispersive_pml', case_dispersive_pml),
         ('3d_mixed',       case_3d_mixed),
         ('steady_state',   case_steady_state)]

for name, case in cases:
    print(f'Testing case: {name}')
    paths = {engine: os.path.join(tempfile.gettempdir(), f'GPU_Engine_{name}_{engine}') for engine in ('basic', 'gpu')}
    logs  = {engine: run_captured(case, paths[engine], engine) for engine in ('basic', 'gpu')}

    assert 'Create FDTD engine (GPU' in logs['gpu'], \
        f'FAIL [{name}]: engine=gpu did not create the GPU engine'
    assert 'Create FDTD engine (GPU' not in logs['basic'], \
        f'FAIL [{name}]: engine=basic created the GPU engine'

    diff, n_probes, n_dumps = compare_outputs(paths['basic'], paths['gpu'])
    print(f'  compared {n_probes} probe files and {n_dumps} field dumps')
    assert not diff, f'FAIL [{name}]: GPU engine differs from the basic engine in: {", ".join(diff)}'
    print('PASS [{}]'.format(name))

print('PASS')
