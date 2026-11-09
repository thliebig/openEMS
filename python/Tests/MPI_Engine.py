# -*- coding: utf-8 -*-
"""
 MPI domain decomposition test — comparison with a serial run

 Each case is written to XML and run by the openEMS binary, once serially and
 once per decomposition with mpirun. Every rank runs the same update on its
 block of the mesh and the halo exchange copies the shared planes exactly, so
 all probe files have to be bit-identical to the serial run, including the
 timestep the end criteria stopped at.

 Together the cases split through PML, Mur ABC, lumped port and RLC element,
 Lorentz/Debye materials and a conducting sheet, in 1D, 2D and 3D
 decompositions with both CPU engines, and stop on the energy end criteria and
 on the steady-state detection. The splits move past the voltage and current
 probes; where no split is possible, all ranks stop with an error. Field dumps
 (TD and FD, HDF5 and VTK, sub-sampled, node and cell interpolated) and a
 NF2FF box cross the splits: rank 0 writes them whole, as a serial run.

 Pass criteria (per case and decomposition)
   the decomposition requested was used
   all probe files, field dumps and far fields bit-identical to the serial run
   no split possible past a probe: exit status != 0 and the error message

 Skipped (with a message) unless the openEMS binary was built with
 -DWITH_MPI=ON and mpirun is found. The binary is $OPENEMS_BIN, or openEMS
 on the PATH.

 Tested with
  - python 3.13
  - openEMS v0.37+, OpenMPI 4.1

"""

import os, re, sys, shutil, subprocess, tempfile
import numpy as np

from CSXCAD import ContinuousStructure
from CSXCAD.CSProperties import CSPropLorentzMaterial, CSPropDebyeMaterial
from openEMS import openEMS
from openEMS.compare import compare_sim_results
from openEMS.ports import LumpedPort

unit = 1e-3   # drawing unit: mm


def find_tools():
    """ the openEMS binary and the mpirun command, None if MPI is not available """
    exe = os.environ.get('OPENEMS_BIN') or shutil.which('openEMS')
    mpirun = shutil.which('mpirun')
    if not exe or not mpirun:
        return None, None, 'openEMS binary or mpirun not found'
    help_text = subprocess.run([exe, '--help'], capture_output=True, text=True).stdout
    if '--mpi-split' not in help_text:
        return None, None, exe + ' was built without MPI (-DWITH_MPI=ON)'
    cmd = [mpirun]
    version = subprocess.run([mpirun, '--version'], capture_output=True, text=True).stdout
    if 'Open MPI' in version or 'OpenRTE' in version:
        cmd += ['--oversubscribe']
        if hasattr(os, 'geteuid') and os.geteuid() == 0:
            cmd += ['--allow-run-as-root']    # CI containers
    return exe, cmd, None


def patch_mesh(CSX, h=1.6):
    mesh = CSX.GetGrid()
    mesh.SetDeltaUnit(unit)
    mesh.AddLine('x', np.concatenate([np.arange(-30, -10, 2.0), np.arange(-10, 10.01, 0.5), np.arange(12, 30.01, 2.0)]))
    mesh.AddLine('y', np.concatenate([np.arange(-24, -8, 2.0), np.arange(-8, 8.01, 0.5), np.arange(10, 24.01, 2.0)]))
    mesh.AddLine('z', np.concatenate([np.arange(-10, 0, 2.0), np.linspace(0, h, 5), np.arange(2.5, 20.01, 1.5)]))


def case_patch():
    """ patch on a lossy substrate, lumped port and series RLC load, PML, energy end criteria """
    h = 1.6
    FDTD = openEMS(NrTS=20000, EndCriteria=1e-4)
    FDTD.SetGaussExcite(4e9, 3.5e9)
    FDTD.SetBoundaryCond(['PML_8'] * 6)
    CSX = ContinuousStructure()
    FDTD.SetCSX(CSX)
    patch_mesh(CSX, h)
    CSX.AddMetal('gnd').AddBox([-20, -16, 0], [20, 16, 0], priority=10)
    CSX.AddMaterial('substrate', epsilon=3.38, kappa=0.002).AddBox([-20, -16, 0], [20, 16, h])
    CSX.AddMetal('patch').AddBox([-8, -6, h], [8, 6, h], priority=10)
    LumpedPort(CSX, 1, 50, [-3, 0, 0], [-3, 0, h], 'z', excite=1)
    rlc = CSX.AddLumpedElement('load', ny='z', caps=True, R=20, C=1e-12, LEtype=0)
    rlc.AddBox([5, 3, 0], [5, 3, h], priority=5)
    CSX.AddProbe('u_load', p_type=0).AddBox([5, 3, 0], [5, 3, h])
    CSX.AddProbe('e_far', p_type=2).AddPoint([20, 15, 10])
    CSX.AddProbe('h_far', p_type=3).AddPoint([-20, -15, 12])
    CSX.AddProbe('e_above', p_type=2).AddPoint([0, 0, 6])
    return FDTD, CSX


def case_materials():
    """ Drude, Lorentz and Debye blocks and a conducting sheet across the splits, Mur ABC """
    FDTD = openEMS(NrTS=1200, EndCriteria=0)
    FDTD.SetGaussExcite(5e9, 4.5e9)
    FDTD.SetBoundaryCond(['MUR'] * 6)
    CSX = ContinuousStructure()
    FDTD.SetCSX(CSX)
    mesh = CSX.GetGrid()
    mesh.SetDeltaUnit(unit)
    for ax in 'xyz':
        mesh.AddLine(ax, np.linspace(-15, 15, 41))
    def lorentz(name, p0, p1, **kw):
        m = CSPropLorentzMaterial(CSX.GetParameterSet(), order=1)
        m.SetName(name)
        m.SetDispersiveMaterialProperty(0, **kw)
        CSX.AddProperty(m)
        m.AddBox(p0, p1, priority=10)
    lorentz('drude', [-8, -8, -8], [3, 3, -2], eps_plasma=5e9, eps_relax=1e-9)
    lorentz('lorentz', [-3, -3, 2], [8, 8, 8], eps_plasma=4e9, eps_pole_freq=3e9, eps_relax=1e-9,
            mue_plasma=3e9, mue_relax=1e-8)
    debye = CSPropDebyeMaterial(CSX.GetParameterSet(), order=2, epsilon=4)
    debye.SetName('debye')
    debye.SetDispersiveMaterialProperty(0, eps_delta=1, eps_relax=4e-11)
    debye.SetDispersiveMaterialProperty(1, eps_delta=0.5, eps_relax=1e-10)
    CSX.AddProperty(debye)
    debye.AddBox([-10, 2, -6], [0, 10, 6], priority=10)
    CSX.AddConductingSheet('sheet', conductivity=1e5, thickness=10e-6).AddBox([-12, -12, 0], [12, 12, 0], priority=20)
    CSX.AddExcitation('dipole', exc_type=0, exc_val=[0, 0, 1]).AddBox([-5.25, 5.25, -3], [-5.25, 5.25, 3])
    for i, p in enumerate([[0.4, 0.4, 0.4], [7, -7, 5], [-11, 11, -11], [11, 11, 11]]):
        CSX.AddProbe('e%d' % i, p_type=2).AddPoint(p)
        CSX.AddProbe('h%d' % i, p_type=3).AddPoint(p)
    return FDTD, CSX


def case_steady_state():
    """ sinusoidal excitation, stopped by the steady-state detection """
    FDTD = openEMS(NrTS=50000, EndCriteria=1e-5)
    FDTD.SetSinusExcite(3e9)
    FDTD.SetBoundaryCond(['PML_8'] * 6)
    CSX = ContinuousStructure()
    FDTD.SetCSX(CSX)
    mesh = CSX.GetGrid()
    mesh.SetDeltaUnit(unit)
    for ax in 'xyz':
        mesh.AddLine(ax, np.linspace(-40, 40, 33))
    CSX.AddMaterial('block', epsilon=4, kappa=0.05).AddBox([-10, -10, -10], [10, 10, 10])
    CSX.AddExcitation('dipole', exc_type=0, exc_val=[0, 0, 1]).AddBox([-20, 0, -5], [-20, 0, 5])
    CSX.AddProbe('e_in', p_type=2).AddPoint([2, 2, 2])
    CSX.AddProbe('e_out', p_type=2).AddPoint([25, -25, 25])
    return FDTD, CSX


def case_probe_long():
    """ a voltage probe from x=-8 to 8 and a current probe below the patch: the splits move past them """
    FDTD, CSX = case_patch()
    CSX.AddProbe('u_long', p_type=0).AddBox([-8, 0, 6], [8, 0, 6])
    CSX.AddProbe('i_patch', p_type=1, norm_dir=2).AddBox([-5, -4, 0.8], [5, 4, 0.8])
    return FDTD, CSX


def case_dumps():
    """ the patch with TD/FD field dumps across the splits and a NF2FF box """
    FDTD, CSX = case_patch()
    CSX.AddDump('Et_vol', dump_type=0, dump_mode=1, file_type=1, sub_sampling=[2, 2, 1]).AddBox([-15, -12, -4], [15, 12, 9])
    CSX.AddDump('Ht_cell', dump_type=1, dump_mode=2, file_type=1, opt_resolution=[1, 1, 1]).AddBox([-20, -16, 0.8], [20, 16, 0.8])
    CSX.AddDump('Jt', dump_type=2, dump_mode=0, file_type=1).AddBox([-20, -16, 0.8], [20, 16, 0.8])
    CSX.AddDump('Ef', dump_type=10, dump_mode=0, file_type=1, frequency=[3e9, 4e9]).AddBox([-10, -8, -2], [10, 8, 4])
    CSX.AddDump('Hf_cell', dump_type=11, dump_mode=2, file_type=1, frequency=[3.5e9]).AddBox([-20, 0, -8], [20, 0, 16])
    CSX.AddDump('Et_vtk', dump_type=0, dump_mode=1, file_type=0).AddBox([-20, -16, 6], [20, 16, 6])
    case_dumps.nf2ff = FDTD.CreateNF2FFBox(start=[-14, -11, -3], stop=[14, 11, 8], frequency=[3.5e9])
    return FDTD, CSX


def post_dumps(sim_path):
    """ the far field from the NF2FF box, into sim_path """
    case_dumps.nf2ff.CalcNF2FF(sim_path, 3.5e9, np.arange(0, 181, 15), [0, 90], read_cached=False)


def case_probe_full():
    """ a voltage probe across the whole mesh between the PMLs in x: no split in x is possible """
    FDTD, CSX = case_patch()
    CSX.AddProbe('u_full', p_type=0).AddBox([-29, 0, 6], [29, 0, 6])
    return FDTD, CSX


def run(exe, mpicmd, xml, sim_path, n=1, split=None, engine=None):
    os.makedirs(sim_path)
    shutil.copy(xml, os.path.join(sim_path, 'sim.xml'))
    cmd = ([] if n == 1 else mpicmd + ['-n', str(n)]) + [exe, 'sim.xml', '--numThreads=1']
    if split:
        cmd.append('--mpi-split=' + split)
    if engine:
        cmd.append('--engine=' + engine)
    p = subprocess.run(cmd, cwd=sim_path, capture_output=True, text=True)
    return p.returncode, p.stdout + p.stderr


def write_xml(case, path):
    FDTD, CSX = case()
    FDTD.Write2XML(path)


exe, mpicmd, why = find_tools()
if exe is None:
    print('MPI_Engine: SKIP, ' + why)
    sys.exit(0)

Sim_Base = os.path.join(tempfile.gettempdir(), 'MPI_Engine')
if os.path.exists(Sim_Base):
    shutil.rmtree(Sim_Base)
os.makedirs(Sim_Base)

#        case                ranks and splits (x,y,z)                     engine
cases = [(case_patch,        [(2, '2,1,1'), (2, '1,2,1'), (4, '2,2,1'), (3, '3,1,1')], None),
         (case_patch,        [(4, '2,2,1')],                                'basic'),
         (case_materials,    [(2, '1,1,2'), (8, '2,2,2'), (4, '0,0,0')],    None),
         (case_materials,    [(6, '1,2,3')],                                'basic'),
         (case_steady_state, [(4, '2,1,2')],                                None),
         (case_probe_long,   [(2, '2,1,1'), (4, '2,2,1')],                  None),
         (case_dumps,        [(2, '2,1,1'), (4, '1,2,2'), (8, '2,2,2')],    None),
         (case_dumps,        [(4, '2,2,1')],                                'basic')]

failures = []
for case, decomps, engine in cases:
    name = case.__name__[5:] + ('_' + engine if engine else '')
    xml = os.path.join(Sim_Base, name + '.xml')
    write_xml(case, xml)
    ref = os.path.join(Sim_Base, name + '_serial')
    rc, out = run(exe, mpicmd, xml, ref, engine=engine)
    assert rc == 0, 'FAIL: serial run of {} failed:\n{}'.format(name, out[-2000:])
    post = globals().get('post_' + case.__name__[5:])
    if post:
        post(ref)
    for n, split in decomps:
        label = '{}_n{}_{}'.format(name, n, split.replace(',', ''))
        path = os.path.join(Sim_Base, label)
        rc, out = run(exe, mpicmd, xml, path, n, split, engine)
        if rc != 0:
            failures.append(label + ': run failed:\n' + out[-2000:])
            continue
        if post:
            post(path)
        m = re.search(r'decomposition into (\d+)x(\d+)x(\d+) blocks', out)
        want = [int(s) for s in split.split(',')]
        if not m or any(w and int(g) != w for g, w in zip(m.groups(), want)):
            failures.append('{}: decomposition not used ({})'.format(label, m.group(0) if m else 'none'))
            continue
        diff = [(f, d) for f, d in compare_sim_results(ref, path) if d != 0]
        print('{:28s} {}  {} probes/datasets {}'.format(label, m.group(0),
              len(compare_sim_results(ref, path)), 'bit-identical' if not diff else 'DIFFER'))
        if diff:
            failures.append('{}: differs from the serial run: {}'.format(label, diff))

# no split in x past a probe across the whole mesh: all ranks stop with an error
xml = os.path.join(Sim_Base, 'probe_full.xml')
write_xml(case_probe_full, xml)
rc, out = run(exe, mpicmd, xml, os.path.join(Sim_Base, 'probe_full_n2'), 2, '2,1,1')
print('{:28s} exit status {}'.format('probe_full_n2_211', rc))
if rc == 0 or 'cannot split direction 0' not in out:
    failures.append('probe_full: expected an error exit, got {}:\n{}'.format(rc, out[-2000:]))

assert not failures, 'FAIL:\n' + '\n'.join(failures)
print('MPI_Engine: all decompositions bit-identical to the serial runs')
shutil.rmtree(Sim_Base)
