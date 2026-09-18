# -*- coding: utf-8 -*-
"""
 Field probe and field dump test — infinitesimal dipole in free space

 A short z-directed soft E-field source radiates in free space (Mur ABC).
 E/H point probes are compared against HDF5 field dumps, time-domain and
 frequency-domain recordings are cross-checked, and the probed near field
 is compared to the analytic Hertzian dipole field.
 (Python port and extension of TESTSUITE/probes/fieldprobes.m)

 Pass criteria
   TD probe vs. TD dump (E, H):  same time axis, max rel. difference < 1e-6
   D-dump = eps0 * E-dump,  B-dump = mue0 * H-dump  (free space, node interpolated),
                                  rel. difference < 1e-5
   FD probe vs. DFT of TD probe:  rel. difference < 1 %
   FD dump  vs. FD probe:         rel. difference < 1e-5
   Symmetry: Ez at (+-r,0,0) and (0,+-r,0) equal,  rel. difference < 1e-5
   Dipole near field:  Ez(r2)/Ez(r1) and Hy(r2)/Hy(r1) within 10 % of the
                       analytic Hertzian dipole ratios (equatorial plane; the
                       Mur ABC is only 20 cells away, inside the reactive near
                       field, which limits the accuracy to a few percent)
   VTK dump files are written

 Tested with
  - python 3.13
  - openEMS v0.0.37+

 (c) 2026 Sean Mollet <sean@malmoset.com>

"""

import os, tempfile, glob
import numpy as np

from CSXCAD  import ContinuousStructure
from openEMS import openEMS
from openEMS.physical_constants import C0, EPS0, MUE0
from openEMS.utilities import DFT_time2freq, HDF5Dump

### Setup the simulation
Sim_Path = os.path.join(tempfile.gettempdir(), 'FieldProbes')

unit     = 1e-3                  # drawing unit: mm
f_max    = 1e9
f_0      = 0.5e9                 # frequency for FD probes / dumps
dip_len  = 6                     # ~ lambda/50 at f_max (mm)
res      = dip_len / 2           # 3 mm mesh
r1       = 5*res                 # probe distances from the dipole (mm)
r2       = 8*res

### FDTD
FDTD = openEMS(NrTS=10000, EndCriteria=1e-6, OverSampling=10)
# band-pass pulse (0.2 .. 1 GHz): a soft source driven with DC content would
# leave a static field behind and spoil the FD evaluation
FDTD.SetGaussExcite(0.6e9, 0.4e9)
FDTD.SetBoundaryCond(['MUR'] * 6)

### CSXCAD geometry & mesh
CSX = ContinuousStructure()
FDTD.SetCSX(CSX)
mesh = CSX.GetGrid()
mesh.SetDeltaUnit(unit)
for ax in 'xyz':
    mesh.AddLine(ax, np.linspace(-20*res, 20*res, 41))

### Infinitesimal dipole: soft Ez source across one cell
dipole = CSX.AddExcitation('infDipole', exc_type=0, exc_val=[0, 0, 1])
dipole.AddBox([0, 0, -dip_len/2], [0, 0, dip_len/2])

### Field dumps on the plane z = 0 (HDF5)
plane_start = [-12*res, -12*res, 0]
plane_stop  = [ 12*res,  12*res, 0]
# no interpolation: raw Yee values, compared against the field probes
for name, dump_type in {'Et': 0, 'Ht': 1}.items():
    dump = CSX.AddDump(name, dump_type=dump_type, dump_mode=0, file_type=1)
    dump.AddBox(plane_start, plane_stop)
# node interpolation: E/D and H/B on a common mesh
for name, dump_type in {'Et_n': 0, 'Dt_n': 4, 'Ht_n': 1, 'Bt_n': 5}.items():
    dump = CSX.AddDump(name, dump_type=dump_type, dump_mode=1, file_type=1)
    dump.AddBox(plane_start, plane_stop)
for name, dump_type in {'Ef': 10, 'Hf': 11}.items():
    dump = CSX.AddDump(name, dump_type=dump_type, dump_mode=0, file_type=1, frequency=[f_0])
    dump.AddBox(plane_start, plane_stop)

# VTK frequency-domain dump (node interpolated) on the same plane
CSX.AddDump('Ef_vtk', dump_type=10, file_type=0, frequency=[f_0]).AddBox(plane_start, plane_stop)

### E/H field probes on the equatorial plane
# E probes sit on mesh nodes, H probes on the dual (Yee) grid half a cell away,
# so that both map unambiguously onto the field dump positions.  The H-dump of
# the z=0 plane is stored on the dual plane z=-res/2.
probe_pos = {'p1': [ r1, 0, 0], 'p2': [-r1, 0, 0],
             'p3': [ 0, r1, 0], 'p4': [ 0, -r1, 0],
             'p5': [ r2, 0, 0]}
h_probe_pos = {name: [pos[0] + res/2, pos[1] + res/2, -res/2] for name, pos in probe_pos.items()}
for name, pos in probe_pos.items():
    CSX.AddProbe('et_' + name, p_type=2).AddPoint(pos)
    CSX.AddProbe('ht_' + name, p_type=3).AddPoint(h_probe_pos[name])
    CSX.AddProbe('ef_' + name, p_type=2, frequency=[f_0]).AddPoint(pos)

### Run the simulation
FDTD.Run(Sim_Path, cleanup=True)

### Post-processing
def load_td_probe(name):
    """ returns (t, field[3, Nt]) of a time-domain field probe """
    data = np.loadtxt(os.path.join(Sim_Path, name), comments='%')
    return data[:, 0], data[:, 1:4].T

def load_fd_probe(name):
    """ returns field[3] (complex) at f_0 of a frequency-domain field probe """
    data = np.atleast_2d(np.loadtxt(os.path.join(Sim_Path, name + '_FD'), comments='%'))
    assert np.isclose(data[0, 0], f_0), f'FAIL: {name}_FD holds f = {data[0, 0]} Hz, expected {f_0} Hz'
    return data[0, 1::2] + 1j*data[0, 2::2]

def dump_at(dump, pos):
    """ index tuple of the dump-mesh node nearest to *pos* (drawing units) """
    return tuple(dump.NearestIndex(n, pos[n]*unit) for n in range(3))

def load_td_dump(name, pos):
    """ returns (t, field[3, Nt]) of a time-domain dump at the node nearest to *pos* """
    with HDF5Dump(os.path.join(Sim_Path, name + '.h5')) as dump:
        idx = dump_at(dump, pos)
        t, val = zip(*[(t, field[(slice(None),) + idx]) for t, field in dump.IterTD()])
    return np.array(t), np.array(val).T

def rel_diff(a, b):
    return np.max(np.abs(a - b)) / np.max(np.abs(b))

### 1. TD probes vs. TD dumps (same Yee positions, dump_mode 0)
for name in probe_pos:
    for probe, dump, pos in (('et_', 'Et', probe_pos[name]), ('ht_', 'Ht', h_probe_pos[name])):
        t_p, val_p = load_td_probe(probe + name)
        t_d, val_d = load_td_dump(dump, pos)
        assert len(t_p) == len(t_d), \
            f'FAIL: {probe}{name} has {len(t_p)} samples, dump {dump} has {len(t_d)}'
        assert np.max(np.abs(t_p - t_d)) < 1e-13, \
            f'FAIL: {probe}{name} and dump {dump} time axes differ'
        err = rel_diff(val_p, val_d)
        assert err < 1e-6, \
            f'FAIL: {probe}{name} differs from dump {dump} by {err:.2e} (rel.), expected < 1e-6'
print('TD probes match TD dumps')

### 2. D = eps0 E and B = mue0 H in free space
with HDF5Dump(os.path.join(Sim_Path, 'Et_n.h5')) as dE, HDF5Dump(os.path.join(Sim_Path, 'Dt_n.h5')) as dD, \
     HDF5Dump(os.path.join(Sim_Path, 'Ht_n.h5')) as dH, HDF5Dump(os.path.join(Sim_Path, 'Bt_n.h5')) as dB:
    t_idx = dE.NumTimesteps // 4   # a timestep with the pulse inside the dump plane
    E = dE.GetFieldAtIndex(t_idx=t_idx); D = dD.GetFieldAtIndex(t_idx=t_idx)
    H = dH.GetFieldAtIndex(t_idx=t_idx); B = dB.GetFieldAtIndex(t_idx=t_idx)
err_D = rel_diff(D, EPS0*E)
err_B = rel_diff(B, MUE0*H)
print(f'max rel. diff. D - eps0 E = {err_D:.2e},  B - mue0 H = {err_B:.2e}')
assert np.max(np.abs(E)) > 0 and np.max(np.abs(H)) > 0, 'FAIL: E/H dump is empty at the chosen timestep'
assert err_D < 1e-5, f'FAIL: D-dump differs from eps0*E by {err_D:.2e} (rel.), expected < 1e-5'
assert err_B < 1e-5, f'FAIL: B-dump differs from mue0*H by {err_B:.2e} (rel.), expected < 1e-5'

### 3. FD probes vs. DFT of TD probes, and FD dump vs. FD probe
with HDF5Dump(os.path.join(Sim_Path, 'Ef.h5')) as dEf:
    Ef = dEf.GetFieldAtFrequency(f_0)
    Ef_idx = {name: dump_at(dEf, pos) for name, pos in probe_pos.items()}

for name in probe_pos:
    t, val = load_td_probe('et_' + name)
    E_dft  = np.array([DFT_time2freq(t, val[n], [f_0])[0] for n in range(3)])
    E_fd   = load_fd_probe('ef_' + name)
    err = rel_diff(E_fd, E_dft)
    assert err < 0.01, \
        f'FAIL: FD probe ef_{name} differs from DFT of TD probe by {err*100:.2f} %, expected < 1 %'
    E_dump = Ef[(slice(None),) + Ef_idx[name]]
    err = rel_diff(E_dump, E_fd)
    assert err < 1e-5, \
        f'FAIL: FD dump Ef differs from FD probe ef_{name} by {err:.2e} (rel.), expected < 1e-5'
print('FD probes match TD probes and FD dump')

### 4. Rotational symmetry of Ez around the dipole axis
Ez = {name: load_fd_probe('ef_' + name)[2] for name in probe_pos}
for name in ('p2', 'p3', 'p4'):
    err = abs(Ez[name] - Ez['p1']) / abs(Ez['p1'])
    assert err < 1e-5, \
        f'FAIL: Ez at {name} differs from Ez at p1 by {err:.2e} (rel.), expected < 1e-5 (symmetry)'
print('Ez is rotationally symmetric')

### 5. Hertzian dipole near field
#   E_theta ~ sin(theta) (1/r^3 + jk/r^2 - k^2/r) exp(-jkr)
#   H_phi   ~ sin(theta) (1/r^2 + jk/r)        exp(-jkr)
# Ez on the equatorial plane is -E_theta; Hy = H_phi cos(phi) at the dual-grid probe position.
k = 2*np.pi*f_0 / C0
def e_z(pos):
    r = np.linalg.norm(pos)*unit
    return -(1/r**3 + 1j*k/r**2 - k**2/r) * np.exp(-1j*k*r)
def h_y(pos):
    x, y, z = np.array(pos)*unit
    r, rho = np.sqrt(x**2 + y**2 + z**2), np.sqrt(x**2 + y**2)
    return (rho/r) * (x/rho) * (1/r**2 + 1j*k/r) * np.exp(-1j*k*r)

def fd_of_td_probe(name, comp):
    t, val = load_td_probe(name)
    return DFT_time2freq(t, val[comp], [f_0])[0]

Hy = {name: fd_of_td_probe('ht_' + name, 1) for name in ('p1', 'p5')}
ratio_E     = Ez['p5'] / Ez['p1']
ratio_E_ana = e_z(probe_pos['p5']) / e_z(probe_pos['p1'])
ratio_H     = Hy['p5'] / Hy['p1']
ratio_H_ana = h_y(h_probe_pos['p5']) / h_y(h_probe_pos['p1'])
err_E = abs(ratio_E / ratio_E_ana - 1)
err_H = abs(ratio_H / ratio_H_ana - 1)
print(f'Ez(r2)/Ez(r1) = {abs(ratio_E):.4f}  (analytic {abs(ratio_E_ana):.4f}),  error {err_E*100:.2f} %')
print(f'Hy(r2)/Hy(r1) = {abs(ratio_H):.4f}  (analytic {abs(ratio_H_ana):.4f}),  error {err_H*100:.2f} %')
assert err_E < 0.10, f'FAIL: Ez distance ratio off by {err_E*100:.1f} % from Hertzian dipole, expected < 10 %'
assert err_H < 0.10, f'FAIL: Hy distance ratio off by {err_H*100:.1f} % from Hertzian dipole, expected < 10 %'

### 6. VTK dump
vtk_files = glob.glob(os.path.join(Sim_Path, 'Ef_vtk*.vtr'))
assert len(vtk_files) > 0, 'FAIL: no VTK files written for the Ef_vtk dump'
assert all(os.path.getsize(f) > 0 for f in vtk_files), 'FAIL: empty VTK dump file'
print(f'{len(vtk_files)} VTK files written')

print('PASS')

if 0:  # set to 1 for debugging plots
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 1, num='Probe vs. dump', tight_layout=True, sharex=True)
    for ax, (probe, dump, comp) in zip(axes, (('et_', 'Et', 2), ('ht_', 'Ht', 1))):
        t_p, val_p = load_td_probe(probe + 'p1')
        t_d, val_d = load_td_dump(dump, probe_pos['p1'])
        ax.plot(t_p*1e9, val_p[comp], 'k-',  linewidth=2, label='probe')
        ax.plot(t_d*1e9, val_d[comp], 'r--', linewidth=1, label='dump')
        ax.set_ylabel(dump)
        ax.grid()
        ax.legend()
    axes[-1].set_xlabel('Time (ns)')
    plt.show()
