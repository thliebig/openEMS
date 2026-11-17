# -*- coding: utf-8 -*-
"""
Material table test — the circuit values from the table and dense must agree bitwise.

The operator keeps the effective material of each node once per distinct value
(the material table) and computes the circuit values from it and the mesh. If
the materials hardly repeat, it falls back to dense circuit values, which normal
models never reach: this test forces each fallback with a low limit
(OPENEMS_EC_TABLE_LIMIT) and compares all probes and field dumps bitwise with
the run on the table.

1. Cartesian, with a lumped R and RC and a conducting sheet, whose circuit
   values go into the table as entries of their own; D-, B- and J-field dumps
   read the material from the table
   a. dense from the start: the limit 1 ends the table during the build
   b. dense after the extensions: one build thread, the limit the size of
      the table before the extensions
2. Cylindrical multi-grid with a D-field dump: the material of the sub-grid's
   region goes into the table, its circuit values stay 0; dense from the start

The field dumps lie away from the lumped elements and the sheet: at those edges
the dumps of the table show the element (the current through a lumped
resistor), the dense ones from the start the material around it.
"""

import os, sys, re, tempfile, ctypes, ctypes.util
import numpy as np

from CSXCAD import ContinuousStructure
from openEMS import openEMS
from openEMS.compare import compare_sim_results


def case_cartesian():
    FDTD = openEMS(NrTS=1500, EndCriteria=0)
    FDTD.SetGaussExcite(10e9, 10e9)
    FDTD.SetBoundaryCond(['PML_8']*6)
    CSX = ContinuousStructure()
    FDTD.SetCSX(CSX)
    # graded along x, the materials average differently at every line
    x = np.unique(np.concatenate([-20 + 40*np.linspace(0, 1, 31)**1.3, [-9, -8, 0, 3, 8, 9]]))
    CSX.DefineGrid({'x': x, 'y': np.arange(-15, 15.1, 1.0), 'z': np.arange(-12, 18.1, 1.0)}, 1e-3)
    CSX.AddMaterial('sub', epsilon=4.3, kappa=0.05).AddBox([-20, -15, -12], [20, 15, -2])
    CSX.AddMaterial('mag', epsilon=2, mue=3, sigma=200).AddBox([-6, -5, -2], [7, 6, 5])
    CSX.AddExcitation('e', 0, [0, 0, 1]).AddBox([0, 0, -2], [0, 0, 2])
    # inside the PML_8 an element of the lumped RLC extension explodes (a separate issue), these are the operator's
    CSX.AddLumpedElement('R', ny=2, caps=True, R=50).AddBox([-9, 4, -2], [-9, 4, 4])
    CSX.AddLumpedElement('RC', ny=2, caps=False, R=100, C=1e-12).AddBox([9, 4, -2], [9, 4, 4])
    CSX.AddConductingSheet('sheet', conductivity=5.8e7, thickness=35e-6).AddBox([-8, -10, 8], [8, 10, 8])
    CSX.AddProbe('ut', p_type=0).AddBox([3, 2, 0], [3, 2, 2])
    CSX.AddProbe('it', p_type=1).AddBox([2, 1, 1], [4, 3, 1])
    for name, t in (('E', 0), ('J', 2), ('D', 4), ('B', 5)):
        CSX.AddDump(name, dump_type=t, file_type=1, dump_mode=0).AddBox([-20, -15, 6], [20, 15, 6])   # between the elements and the sheet
    return FDTD


def case_multigrid():
    FDTD = openEMS(CoordSystem=1, NrTS=1500, EndCriteria=0)
    FDTD.SetMultiGrid([18])
    FDTD.SetGaussExcite(3e9, 3e9)
    FDTD.SetBoundaryCond(['PEC', 'PEC', 'PEC', 'PEC', 'PML_8', 'PML_8'])
    CSX = ContinuousStructure(CoordSystem=1)
    FDTD.SetCSX(CSX)
    mesh = CSX.GetGrid()
    mesh.SetDeltaUnit(1e-3)
    mesh.AddLine('r', np.linspace(0, 40, 41))
    mesh.AddLine('a', np.linspace(-np.pi, np.pi, 33))
    mesh.AddLine('z', np.linspace(-20, 20, 41))
    CSX.AddMaterial('ring', epsilon=3, kappa=0.02).AddCylindricalShell([0, 0, -5], [0, 0, 5], 20, 8)
    CSX.AddMaterial('core', epsilon=2.2, mue=1.5).AddCylinder([0, 0, -10], [0, 0, 10], 6)
    CSX.AddExcitation('e', 0, [0, 0, 1]).AddBox([12, 0, -2], [12, 0, 2])
    CSX.AddProbe('ut', p_type=0).AddBox([15, 0.5, -2], [15, 0.5, 2])
    CSX.AddDump('D', dump_type=4, file_type=1, dump_mode=0).AddBox([0, -np.pi, 1], [40, np.pi, 1])
    return FDTD


def run(case, label, limit=None, threads=0):
    """ Run a case, with the material table limited to \a limit entries; return its path and the console output """
    Sim_Path = os.path.join(tempfile.gettempdir(), 'Material_Table_' + label)
    if limit is None:
        os.environ.pop('OPENEMS_EC_TABLE_LIMIT', None)
    else:
        os.environ['OPENEMS_EC_TABLE_LIMIT'] = str(limit)
    log_file = Sim_Path + '.log'
    sys.stdout.flush()
    saved = os.dup(1)
    try:
        with open(log_file, 'w') as log:
            os.dup2(log.fileno(), 1)
            try:
                case().Run(Sim_Path, cleanup=True, verbose=1, numThreads=threads)
            finally:
                try:   # flush the C stdio buffer of the openEMS library
                    ctypes.CDLL(ctypes.util.find_library('c')).fflush(None)
                except (OSError, AttributeError, TypeError):
                    pass
                os.dup2(saved, 1)
                os.close(saved)
    finally:
        os.environ.pop('OPENEMS_EC_TABLE_LIMIT', None)
    with open(log_file) as log:
        text = log.read()
    # all timesteps: a NaN energy ends a run without end criteria early, and silently
    m = re.search(r'Time for (\d+) iterations', text)
    assert m and int(m.group(1))==1500, f'FAIL [{label}]: the run ended early: ' + (m.group(0) if m else 'no timing line')
    return Sim_Path, text


def check_same(path_ref, path, label):
    files = [f for f in os.listdir(path_ref) if os.path.isfile(os.path.join(path_ref, f))]
    assert any(np.any(np.loadtxt(os.path.join(path_ref, f), comments='%')[:, 1:]) for f in files if '.' not in f), \
        f'FAIL [{label}]: the probes are zero, the comparison would test nothing'
    diff = compare_sim_results(path_ref, path)
    bad = [(n, d) for n, d in diff if d != 0]
    print(f'  {label}: {len(diff)} probes and dump datasets, {"all bitwise equal" if not bad else str(len(bad)) + " differ"}')
    assert not bad, f'FAIL [{label}]: differs from the run on the material table: {bad[:5]}{" ..." if len(bad)>5 else ""}'


# 1. Cartesian
path_table, log = run(case_cartesian, 'cart_table', threads=1)
m = re.search(r'Operator: (\d+) distinct materials', log)
assert m and 'go dense' not in log, 'FAIL: the reference run is not on the material table'
n_table = int(m.group(1))
m = re.search(r'Operator: (\d+) entries in the material table with the extensions', log)
assert m and int(m.group(1)) > n_table, 'FAIL: the lumped elements and the sheet add no entries to the table'
print(f'material table: {n_table} entries, {m.group(1)} with the extensions')

path, log = run(case_cartesian, 'cart_dense_build', limit=1)
assert 'go dense' in log, 'FAIL [dense from the start]: the run stayed on the material table'
check_same(path_table, path, 'dense from the start')

path, log = run(case_cartesian, 'cart_dense_ext', limit=n_table, threads=1)
assert re.search(r'Operator: \d+ distinct materials, more than the material table is worth', log), \
    'FAIL [dense after the extensions]: the run did not go dense when the extensions added entries'
check_same(path_table, path, 'dense after the extensions')

# 2. cylindrical multi-grid
path_table, log = run(case_multigrid, 'mg_table')
assert 'go dense' not in log, 'FAIL: the multi-grid reference run is not on the material table'
path, log = run(case_multigrid, 'mg_dense', limit=1)
assert 'go dense' in log, 'FAIL [multi-grid dense]: the run stayed on the material table'
check_same(path_table, path, 'multi-grid, dense')

print('PASS')
