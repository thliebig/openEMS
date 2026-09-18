# -*- coding: utf-8 -*-
"""
 Cylindrical coordinates test — PEC pillbox cavity resonator

 A closed cylindrical cavity (radius R, height h) is simulated on a full
 (closed) cylindrical mesh including the r=0 axis.  An off-axis, z-directed
 line source spanning the full height excites only TM_mn0 modes, with

     f_mn0 = x_mn * c0 / (2 pi R)      (x_mn: n-th zero of J_m)

 so below the first TE/TM_mn1 mode only TM010 and TM110 resonate in the
 evaluated band.  The test is run on a plain cylindrical mesh and repeated
 with a cylindrical multi-grid (coarser alpha resolution near the axis).

 Pass criteria (per mesh)
   TM010 and TM110 resonance frequencies within 1 % of the analytic values
   TM010 radial field profile (alpha averaged FD field dump) matches
     J0(x_01 r / R) within 0.03 (normalised)
   plain and multi-grid resonance frequencies agree within 0.5 %

 Tested with
  - python 3.13
  - openEMS v0.0.37+

 (c) 2026 Sean Mollet <sean@malmoset.com>

"""

import os, tempfile
import numpy as np

from CSXCAD  import ContinuousStructure
from openEMS import openEMS
from openEMS.physical_constants import C0
from openEMS.utilities import DFT_time2freq, HDF5Dump

unit   = 1e-3     # drawing unit: mm
R      = 50       # cavity radius (mm)
height = 30       # cavity height (mm), TM/TE_mn1 modes lie above 5 GHz
dr     = 2.5      # radial mesh resolution (mm)
dz     = 2.5

x_01 = 2.404826   # first zero of J0
x_11 = 3.831706   # first zero of J1
f_010 = x_01 * C0 / (2*np.pi*R*unit)   # 2.295 GHz
f_110 = x_11 * C0 / (2*np.pi*R*unit)   # 3.657 GHz

f_start = 1.5e9
f_stop  = 4.5e9
freq    = np.linspace(f_start, f_stop, 3001)

r_src   = 7*dr      # source position (r, alpha=0), on a mesh line
r_probe = 10*dr     # probe position  (r, alpha=pi/4)
a_probe = np.pi/4


def bessel_j0(x):
    """ J0(x) = 1/pi * int_0^pi cos(x sin(t)) dt """
    t = np.linspace(0, np.pi, 2001)
    y = np.cos(np.outer(x, np.sin(t)))
    return np.sum(0.5*(y[:, 1:] + y[:, :-1]), axis=1) * (t[1] - t[0]) / np.pi


def peak_near(f, spec, f_ana):
    """ frequency of the spectral maximum within +/- 5 % of f_ana """
    mask = np.abs(f - f_ana) < 0.05*f_ana
    return f[mask][np.argmax(spec[mask])]


def run(label, n_alpha, multi_grid=None):
    Sim_Path = os.path.join(tempfile.gettempdir(), 'CylCavity_' + label)

    FDTD = openEMS(CoordSystem=1, NrTS=40000, EndCriteria=0)
    if multi_grid is not None:
        FDTD.SetMultiGrid(multi_grid)
    FDTD.SetGaussExcite(0.5*(f_start + f_stop), 0.5*(f_stop - f_start))
    FDTD.SetBoundaryCond(['PEC'] * 6)   # r=0 and the closed alpha direction need no BC

    CSX = ContinuousStructure(CoordSystem=1)
    FDTD.SetCSX(CSX)
    mesh = CSX.GetGrid()
    mesh.SetDeltaUnit(unit)
    mesh.AddLine('r', np.arange(0, R + dr/2, dr))
    mesh.AddLine('a', (np.arange(n_alpha + 1) - n_alpha/2) * 2*np.pi/n_alpha)   # exact line at alpha=0
    mesh.AddLine('z', np.arange(0, height + dz/2, dz))

    src = CSX.AddExcitation('line_source', exc_type=0, exc_val=[0, 0, 1])
    src.AddBox([r_src, 0, 0], [r_src, 0, height])

    CSX.AddProbe('et', p_type=2).AddPoint([r_probe, a_probe, height/2])

    dump = CSX.AddDump('Ef', dump_type=10, dump_mode=1, file_type=1, frequency=[f_010])
    dump.AddBox([0, -np.pi, height/2], [R, np.pi, height/2])

    FDTD.Run(Sim_Path, cleanup=True)

    # probe columns: t/s, E_rho, E_alpha, E_z
    data = np.loadtxt(os.path.join(Sim_Path, 'et'), comments='%')
    spec = np.abs(DFT_time2freq(data[:, 0], data[:, 3], freq))

    f_res = {'TM010': peak_near(freq, spec, f_010),
             'TM110': peak_near(freq, spec, f_110)}
    for mode, f_ana in (('TM010', f_010), ('TM110', f_110)):
        err = f_res[mode]/f_ana - 1
        print(f'  {mode}: {f_res[mode]/1e9:.4f} GHz  (analytic {f_ana/1e9:.4f} GHz),  error {err*100:+.2f} %')
        assert abs(err) < 0.01, \
            f'FAIL [{label}]: {mode} at {f_res[mode]/1e9:.4f} GHz, expected {f_ana/1e9:.4f} GHz (+/- 1 %)'

    # TM010 radial profile: average Ez over alpha, compare to J0
    with HDF5Dump(os.path.join(Sim_Path, 'Ef.h5')) as fd:
        r_lines = fd.GetMesh()['lines'][0] / unit
        Ez = fd.GetFieldAtFrequency(f_010, component=2)[:, :, 0]   # (r, alpha)
    prof     = np.abs(np.mean(Ez, axis=1))
    prof    /= prof[0]
    prof_ana = bessel_j0(x_01 * r_lines / R)
    prof_err = np.max(np.abs(prof - prof_ana))
    print(f'  TM010 radial profile: max deviation from J0 = {prof_err:.4f}')
    assert prof_err < 0.03, \
        f'FAIL [{label}]: TM010 radial profile deviates {prof_err:.4f} from J0, expected < 0.03'

    return f_res, spec


print('Plain cylindrical mesh')
f_plain, spec_plain = run('plain', n_alpha=24)
print('PASS [plain]')

print('Cylindrical multi-grid mesh (alpha resolution halved inside r = R/2)')
f_mg, spec_mg = run('multigrid', n_alpha=48, multi_grid=[R/2])
print('PASS [multigrid]')

for mode in f_plain:
    diff = abs(f_mg[mode]/f_plain[mode] - 1)
    assert diff < 0.005, \
        f'FAIL: {mode} differs by {diff*100:.2f} % between plain and multi-grid mesh, expected < 0.5 %'

print('PASS')

if 0:  # set to 1 for debugging plots
    import matplotlib.pyplot as plt

    fig, axis = plt.subplots(num='Cavity spectrum', tight_layout=True)
    axis.semilogy(freq/1e9, spec_plain, 'k-',  linewidth=2, label='plain')
    axis.semilogy(freq/1e9, spec_mg,    'r--', linewidth=2, label='multi-grid')
    for f_ana in (f_010, f_110):
        axis.axvline(f_ana/1e9, color='g', linestyle=':')
    axis.grid()
    axis.set_xmargin(0)
    axis.set_xlabel('Frequency (GHz)')
    axis.set_ylabel('|Ez| (a.u.)')
    axis.legend()
    plt.show()
