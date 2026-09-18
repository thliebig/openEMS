# -*- coding: utf-8 -*-
"""
 GPU engine test — comparison with the basic engine

 Each case runs with engine='basic', engine='gpu-reference' and engine='gpu'
 and compares all probe files and all HDF5 field dumps.

 The reference backend performs the main FDTD updates on the CPU on its own
 copy of the fields, exactly like the basic engine, so its results must be
 bit-identical. A device backend (engine='gpu': Metal on macOS, CUDA) runs on
 the GPU, whose float arithmetic may differ in rounding (e.g. fused
 multiply-add), so it is compared with a tolerance.

 Together the cases use the excitation, UPML, Mur ABC, Lorentz material,
 lumped RLC, conducting sheet, TF/SF, local absorber, steady-state and
 cylinder extensions, and cylindrical meshes with one and two multi-grid
 levels. For cylindrical meshes the engine choice 'basic' has no effect, the
 CPU reference is the cylindrical (multithreaded) engine.

 Pass criteria (per case)
   the requested backend was created (not a silent fallback)
   device backend: all extensions run on the device (backends in FULL_DEVICE_BACKENDS)
   reference backend: all probe data and field dumps bit-identical
   device backend: max. deviation < 1e-4 of the peak value, per probe file and per field dump
     (skipped without a GPU device)

 Tested with
  - python 3.13
  - openEMS v0.0.37+

 (c) 2026 Sean Mollet <sean@malmoset.com>

"""

import os, re, sys, glob, tempfile, ctypes, ctypes.util
import numpy as np
import h5py

from CSXCAD  import ContinuousStructure
from CSXCAD.CSProperties import CSPropLorentzMaterial, CSPropDebyeMaterial
from CSXCAD.CSProperties import ABCtype
from openEMS.physical_constants import C0
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


def case_excitation():
    """ PEC cavity: only the excitation extension, with overlapping sources """
    FDTD = openEMS(NrTS=600, EndCriteria=0)
    FDTD.SetGaussExcite(5.5e9, 4.5e9)
    FDTD.SetBoundaryCond(['PEC'] * 6)
    CSX = ContinuousStructure()
    FDTD.SetCSX(CSX)
    mesh = CSX.GetGrid()
    mesh.SetDeltaUnit(unit)
    for ax in 'xyz':
        mesh.AddLine(ax, np.arange(-10, 10.5, 1))
    CSX.AddExcitation('e_soft', exc_type=0, exc_val=[1, 0, 0]).AddBox([0, -2, 2], [4, 2, 2])
    CSX.AddExcitation('e_soft2', exc_type=0, exc_val=[0.5, 0, 0], delay=0.1e-9).AddBox([2, 0, 2], [6, 0, 2])  # shares edges
    CSX.AddExcitation('h_soft', exc_type=2, exc_val=[0, 0, 1]).AddBox([2, 2, -4], [4, 4, -2])
    CSX.AddProbe('et', p_type=2).AddPoint([5, 5, 5])
    CSX.AddProbe('ht', p_type=3).AddPoint([-5, 5, -5])
    CSX.AddDump('Et', dump_type=0, file_type=1).AddBox([-10, -10, 0], [10, 10, 0])
    return FDTD, CSX


def case_pml():
    """ free space with PML on all sides: excitation and UPML extensions """
    FDTD = openEMS(NrTS=700, EndCriteria=0)
    FDTD.SetGaussExcite(5.5e9, 4.5e9)
    FDTD.SetBoundaryCond(['PML_8'] * 6)
    CSX = ContinuousStructure()
    FDTD.SetCSX(CSX)
    mesh = CSX.GetGrid()
    mesh.SetDeltaUnit(unit)
    for ax in 'xyz':
        mesh.AddLine(ax, np.arange(-15, 15.5, 1))
    CSX.AddExcitation('dipole', exc_type=0, exc_val=[0, 0, 1]).AddBox([0, 0, -1], [0, 0, 1])
    CSX.AddProbe('et', p_type=2).AddPoint([4, 3, 2])
    CSX.AddProbe('ht', p_type=3).AddPoint([-3, 5, 0])
    CSX.AddDump('Et', dump_type=0, file_type=1).AddBox([-15, -15, 0], [15, 15, 0])
    return FDTD, CSX


def case_mur():
    """ free space with Mur ABC on all sides and a source on a boundary plane: excitation and Mur extensions """
    FDTD = openEMS(NrTS=700, EndCriteria=0)
    FDTD.SetGaussExcite(5.5e9, 4.5e9)
    FDTD.SetBoundaryCond(['MUR'] * 6)
    CSX = ContinuousStructure()
    FDTD.SetCSX(CSX)
    mesh = CSX.GetGrid()
    mesh.SetDeltaUnit(unit)
    for ax in 'xyz':
        mesh.AddLine(ax, np.arange(-10, 10.5, 1))
    CSX.AddExcitation('dipole', exc_type=0, exc_val=[0, 0, 1]).AddBox([0, 0, -1], [0, 0, 1])
    # a source on the x-max plane delays that Mur ABC until the excitation is done
    CSX.AddExcitation('wall', exc_type=0, exc_val=[0, 1, 0]).AddBox([10, -2, 0], [10, 2, 0])
    CSX.AddProbe('et', p_type=2).AddPoint([4, 3, 2])
    CSX.AddProbe('ht', p_type=3).AddPoint([-3, 5, 0])
    CSX.AddDump('Et', dump_type=0, file_type=1).AddBox([-10, -10, 0], [10, 10, 0])
    return FDTD, CSX


def case_materials():
    """ Drude, Lorentz and Debye materials, a magnetic Drude material and a conducting sheet in a PML channel """
    FDTD = openEMS(NrTS=3000, EndCriteria=0)
    CSX = ContinuousStructure()
    FDTD.SetCSX(CSX)
    CSX.GetGrid().SetDeltaUnit(unit)
    channel_1d(FDTD, CSX)
    def lorentz(name, z0, z1, **kw):
        m = CSPropLorentzMaterial(CSX.GetParameterSet(), order=1)
        m.SetName(name)
        m.SetDispersiveMaterialProperty(0, **kw)
        CSX.AddProperty(m)
        m.AddBox([0, 0, z0], [1, 1, z1], priority=10)
    lorentz('drude', 25, 30, eps_plasma=5e9, eps_relax=1e-9)
    lorentz('lorentz', 32, 37, eps_plasma=4e9, eps_pole_freq=3e9, eps_relax=1e-9)
    lorentz('double_drude', 39, 44, eps_plasma=5e9, eps_relax=1e-8, mue_plasma=5e9, mue_relax=1e-8)
    debye = CSPropDebyeMaterial(CSX.GetParameterSet(), order=1, epsilon=4)
    debye.SetName('debye')
    debye.SetDispersiveMaterialProperty(0, eps_delta=1, eps_relax=4e-11)
    CSX.AddProperty(debye)
    debye.AddBox([0, 0, 46], [1, 1, 51], priority=10)
    CSX.AddConductingSheet('sheet', conductivity=1e5, thickness=10e-6).AddBox([0, 0, 55], [1, 1, 55], priority=10)
    return FDTD, CSX


def case_lumped():
    """ lumped port with series and parallel RLC elements in a Mur box: excitation, lumped RLC and Mur extensions """
    FDTD = openEMS(NrTS=1500, EndCriteria=0)
    FDTD.SetGaussExcite(5.5e9, 4.5e9)
    FDTD.SetBoundaryCond(['MUR'] * 6)
    CSX = ContinuousStructure()
    FDTD.SetCSX(CSX)
    mesh = CSX.GetGrid()
    mesh.SetDeltaUnit(unit)
    for ax in 'xyz':
        mesh.AddLine(ax, np.arange(-10, 10.5, 1))
    LumpedPort(CSX, 1, 50, [-4, 0, 0], [-4, 0, 2], 'z', excite=1)
    ser = CSX.AddLumpedElement('ser_rlc', ny='z', caps=False, R=10, L=1e-9, C=1e-12, LEtype=1)
    ser.AddBox([0, 0, 0], [0, 0, 2], priority=10)
    par = CSX.AddLumpedElement('par_rlc', ny='z', caps=False, R=200, L=2e-9, C=0.5e-12, LEtype=0)
    par.AddBox([4, 0, 0], [4, 0, 2], priority=10)
    wire = CSX.AddMetal('wire')
    wire.AddCurve([[-4, 4], [0, 0], [0, 0]])
    wire.AddCurve([[-4, 4], [0, 0], [2, 2]])
    CSX.AddProbe('et', p_type=2).AddPoint([0, 5, 0])
    return FDTD, CSX


def case_tfsf():
    """ oblique plane wave on a PEC sphere in a PML box: excitation, TF/SF and UPML extensions """
    FDTD = openEMS(NrTS=700, EndCriteria=0)
    FDTD.SetGaussExcite(5.5e9, 4.5e9)
    FDTD.SetBoundaryCond(['PML_8'] * 6)
    CSX = ContinuousStructure()
    FDTD.SetCSX(CSX)
    mesh = CSX.GetGrid()
    mesh.SetDeltaUnit(unit)
    for ax in 'xyz':
        mesh.AddLine(ax, np.arange(-15, 15.5, 1))
    k_dir = np.array([1, 2, 3]) / np.sqrt(14)
    pw = CSX.AddExcitation('plane_wave', exc_type=10, exc_val=[2, -1, 0])
    pw.SetPropagationDir(k_dir)
    pw.SetFrequency(5e9)
    pw.AddBox([-6, -6, -6], [6, 6, 6])
    CSX.AddMetal('sphere').AddSphere(priority=10, center=[0, 0, 0], radius=3)
    CSX.AddProbe('et_in', p_type=2).AddPoint([4, -3, 2])
    CSX.AddProbe('et_out', p_type=2).AddPoint([-10, 1, 3])
    CSX.AddDump('Et', dump_type=0, file_type=1).AddBox([-15, -15, 0], [15, 15, 0])
    return FDTD, CSX


def case_absorbers():
    """ PEC-terminated channel with local absorbing sheets (Mur and Mur with super-absorption) """
    FDTD = openEMS(NrTS=3000, EndCriteria=0)
    FDTD.SetGaussExcite(5.5e9, 4.5e9)
    FDTD.SetBoundaryCond(['PEC', 'PEC', 'PMC', 'PMC', 'PEC', 'PEC'])
    CSX = ContinuousStructure()
    FDTD.SetCSX(CSX)
    mesh = CSX.GetGrid()
    mesh.SetDeltaUnit(unit)
    mesh.AddLine('x', [0, 0.5, 1])
    mesh.AddLine('y', [0, 0.5, 1])
    mesh.AddLine('z', np.arange(0, 100.5, 0.5))
    CSX.AddExcitation('plane', exc_type=0, exc_val=[1, 0, 0]).AddBox([0, 0, 40], [1, 1, 40])
    CSX.AddAbsorbingBC('abs_low', NormalSignPositive=False, AbsorbingBoundaryType=ABCtype.MUR_1ST,
                       PhaseVelocity=C0).AddBox([0, 0, 5], [1, 1, 5], priority=6)
    CSX.AddAbsorbingBC('abs_high', NormalSignPositive=True, AbsorbingBoundaryType=ABCtype.MUR_1ST_SA,
                       PhaseVelocity=C0).AddBox([0, 0, 95], [1, 1, 95], priority=6)
    CSX.AddProbe('et', p_type=2).AddPoint([0.5, 0.5, 70])
    CSX.AddProbe('ht', p_type=3).AddPoint([0.5, 0.5, 20])
    return FDTD, CSX


def cylinder_mesh(FDTD, alpha, r0, r1, z_pad=0):
    """ cylindrical mesh, the CPU reference is the cylindrical engine (engine='basic' has no effect);
        z_pad extends z beyond 0..30 on both sides, e.g. to make room for a PML """
    CSX = ContinuousStructure(CoordSystem=1)
    FDTD.SetCSX(CSX)
    mesh = CSX.GetGrid()
    mesh.SetDeltaUnit(unit)
    mesh.AddLine('r', np.arange(r0, r1 + 1, 2))
    mesh.AddLine('a', alpha)
    mesh.AddLine('z', np.arange(-z_pad, 30 + z_pad + 0.5, 2))
    return CSX


def case_cylinder_closed():
    """ closed cylindrical mesh including r=0: excitation and cylinder extensions """
    FDTD = openEMS(CoordSystem=1, NrTS=1500, EndCriteria=0)
    FDTD.SetGaussExcite(3e9, 2e9)
    FDTD.SetBoundaryCond(['PEC'] * 6)
    CSX = cylinder_mesh(FDTD, (np.arange(25) - 12) * 2*np.pi/24, 0, 40)
    CSX.AddExcitation('line', exc_type=0, exc_val=[0, 0, 1]).AddBox([14, 0, 0], [14, 0, 30])
    CSX.AddExcitation('radial', exc_type=0, exc_val=[1, 0, 0]).AddBox([6, np.pi/2, 10], [10, np.pi/2, 10])
    CSX.AddProbe('et_axis', p_type=2).AddPoint([0, 0, 16])
    CSX.AddProbe('et', p_type=2).AddPoint([20, np.pi/4, 14])
    CSX.AddProbe('ht', p_type=3).AddPoint([10, -np.pi/3, 8])
    CSX.AddDump('Et', dump_type=0, file_type=1).AddBox([0, -np.pi, 14], [40, np.pi, 14])
    return FDTD, CSX


def case_cylinder_wedge():
    """ open alpha wedge with r>0 and PML in z: excitation, UPML and (inactive) cylinder extensions """
    FDTD = openEMS(CoordSystem=1, NrTS=1000, EndCriteria=0)
    FDTD.SetGaussExcite(3e9, 2e9)
    FDTD.SetBoundaryCond(['PEC', 'PEC', 'PEC', 'PEC', 'PML_8', 'PML_8'])
    CSX = cylinder_mesh(FDTD, np.linspace(-np.pi/4, np.pi/4, 13), 10, 40, z_pad=16)
    CSX.AddExcitation('coax', exc_type=0, exc_val=[1, 0, 0]).AddBox([10, -np.pi/4, 12], [40, np.pi/4, 12])
    CSX.AddProbe('et', p_type=2).AddPoint([20, 0, 20])
    CSX.AddProbe('ht', p_type=3).AddPoint([30, np.pi/8, 6])
    return FDTD, CSX


def multigrid(radii, alpha, r1=40, z_bc='PEC'):
    """ cylindrical multi-grid with sources and probes inside and outside the sub-grids """
    FDTD = openEMS(CoordSystem=1, NrTS=1200, EndCriteria=0, MultiGrid=radii)
    FDTD.SetGaussExcite(3e9, 2e9)
    FDTD.SetBoundaryCond(['PEC', 'PEC', 'PEC', 'PEC', z_bc, z_bc])
    CSX = cylinder_mesh(FDTD, alpha, 0, r1, z_pad=16 if z_bc.startswith('PML') else 0)
    CSX.AddExcitation('inner', exc_type=0, exc_val=[0, 0, 1]).AddBox([4, 0, 0], [4, 0, 30])
    CSX.AddExcitation('outer', exc_type=0, exc_val=[1, 0, 0]).AddBox([26, alpha[3], 12], [32, alpha[3], 12])
    CSX.AddProbe('et_axis', p_type=2).AddPoint([0, 0, 16])
    CSX.AddProbe('et_inner', p_type=2).AddPoint([6, alpha[len(alpha)//3], 14])
    CSX.AddProbe('ht_inner', p_type=3).AddPoint([8, alpha[len(alpha)//4], 8])
    CSX.AddProbe('et_outer', p_type=2).AddPoint([30, alpha[2*len(alpha)//3], 14])
    CSX.AddDump('Et', dump_type=0, file_type=1).AddBox([0, alpha[0], 14], [r1, alpha[-1], 14])
    return FDTD, CSX


def case_multigrid():
    """ closed cylindrical mesh with one multi-grid level """
    return multigrid([14], (np.arange(49) - 24) * 2*np.pi/48)


def case_multigrid2():
    """ closed cylindrical mesh with two nested multi-grid levels """
    return multigrid([10, 20], (np.arange(49) - 24) * 2*np.pi/48)


def case_multigrid_wedge():
    """ open alpha wedge with one multi-grid level and PML in z """
    return multigrid([14], np.linspace(-np.pi/2, np.pi/2, 25), z_bc='PML_8')


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


def deviation(a, b, peak=None):
    """ max. deviation of b from a, relative to peak (default: the peak of a), 0 if identical """
    a = np.asarray(a); b = np.asarray(b)
    if a.shape != b.shape:
        return np.inf
    if np.array_equal(a, b):
        return 0.0
    if peak is None:
        peak = np.max(np.abs(a))
    return np.max(np.abs(a - b)) / (peak if peak > 0 else 1.0)


def compare_outputs(path_a, path_b, rtol=0):
    """ compare all probe files and HDF5 dumps, return a list of (name, deviation) above rtol """
    diff = []
    files = [f for f in os.listdir(path_a) if os.path.isfile(os.path.join(path_a, f))]
    probes = [f for f in files if '.' not in f]   # probe and port files have no extension
    dumps  = [f for f in files if f.endswith('.h5')]
    assert probes, f'FAIL: no probe files in {path_a}'
    for f in probes:
        a = np.loadtxt(os.path.join(path_a, f), comments='%')
        b = np.loadtxt(os.path.join(path_b, f), comments='%')
        diff.append((f, deviation(a, b)))
    for f in dumps:
        with h5py.File(os.path.join(path_a, f), 'r') as a, h5py.File(os.path.join(path_b, f), 'r') as b:
            names = []
            a.visititems(lambda name, obj: names.append(name) if isinstance(obj, h5py.Dataset) else None)
            # field data relative to the peak of the whole dump, not of each (possibly decayed) timestep
            fields = [n for n in names if n.startswith('FieldData')]
            peak = max([np.max(np.abs(a[n][()])) for n in fields] or [0])
            for n in names:
                if n not in b:
                    diff.append((f'{f}:{n}', np.inf))
                else:
                    diff.append((f'{f}:{n}', deviation(a[n][()], b[n][()], peak if n in fields else None)))
    worst = max(d for _, d in diff)
    return [(n, d) for n, d in diff if d > rtol], len(probes), len(dumps), worst


# (name, case, all extensions have a Metal implementation)
cases = [('excitation',     case_excitation,     True),
         ('pml',            case_pml,            True),
         ('mur',            case_mur,            True),
         ('materials',      case_materials,      True),
         ('lumped',         case_lumped,         True),
         ('tfsf',           case_tfsf,           True),
         ('absorbers',      case_absorbers,      True),
         ('cylinder_closed', case_cylinder_closed, True),
         ('cylinder_wedge', case_cylinder_wedge,  True),
         ('multigrid',      case_multigrid,       True),
         ('multigrid2',     case_multigrid2,      True),
         ('multigrid_wedge', case_multigrid_wedge, True),
         ('dispersive_pml', case_dispersive_pml, True),
         ('3d_mixed',       case_3d_mixed,       True),
         ('steady_state',   case_steady_state,   True)]

DEVICE_RTOL = 1e-4
# device backends with a device implementation of every extension
FULL_DEVICE_BACKENDS = ('Metal',)
engines = ('basic', 'gpu-reference', 'gpu')

for name, case, on_device in cases:
    print(f'Testing case: {name}')
    paths = {engine: os.path.join(tempfile.gettempdir(), f'GPU_Engine_{name}_{engine}') for engine in engines}
    logs  = {engine: run_captured(case, paths[engine], engine) for engine in engines}

    assert 'Create FDTD engine (GPU' not in logs['basic'], \
        f'FAIL [{name}]: engine=basic created the GPU engine'
    assert 'Create FDTD engine (GPU, backend: reference' in logs['gpu-reference'], \
        f'FAIL [{name}]: engine=gpu-reference did not create the GPU engine with the reference backend'

    diff, n_probes, n_dumps, _ = compare_outputs(paths['basic'], paths['gpu-reference'])
    print(f'  reference backend: compared {n_probes} probe files and {n_dumps} field dumps')
    assert not diff, f'FAIL [{name}]: reference backend differs from the basic engine in: {", ".join(n for n, _ in diff)}'

    backend = re.search(r'Create FDTD engine \(GPU, backend: (\w+)', logs['gpu'])
    backend = backend.group(1) if backend else None
    if backend and backend!='reference':
        if on_device and backend in FULL_DEVICE_BACKENDS:
            assert ('Engine_GPU: all extensions run on the device' in logs['gpu']) and ('host fallback' not in logs['gpu']), \
                f'FAIL [{name}]: the {backend} backend did not run all extensions on the device'
        diff, _, _, worst = compare_outputs(paths['basic'], paths['gpu'], rtol=DEVICE_RTOL)
        print(f'  {backend} backend: max. deviation {worst:.1e} of the peak value')
        assert not diff, f'FAIL [{name}]: {backend} backend deviates from the basic engine: ' + \
            ', '.join(f'{n} ({d:.1e})' for n, d in diff)
    else:
        print('  no GPU device, device backend not tested')
    print('PASS [{}]'.format(name))

print('PASS')
