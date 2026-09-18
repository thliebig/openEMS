# -*- coding: utf-8 -*-
"""
 Near-field to far-field (NF2FF) test — half-wave dipole and quarter-wave monopole

 1. Half-wave dipole in free space (PML on all sides), broadband (time-domain)
    NF2FF recording box.  The far field is compared to the analytic thin
    half-wave dipole:

        E(theta) ~ cos(pi/2 cos(theta)) / sin(theta),   D = 1.64 (2.15 dBi)

 2. Quarter-wave monopole on an infinite ground plane (PEC boundary at z=0),
    single-frequency (frequency-domain) NF2FF recording box with PEC mirroring
    of the missing box face.  Image theory: same pattern in the upper half
    space, twice the directivity of the dipole, D = 3.28 (5.16 dBi).

 The radiated power is also obtained by integrating the far-field radiation
 intensity over the (half) sphere.  For the monopole, directivity and
 efficiency are taken from this integral over the real upper half space,
 since the nf2ff P_rad / Dmax include the mirrored (image) half space.

 Pass criteria (at the S11 minimum)
   Dmax within 5 % of the analytic value (1.64 / 3.28)
   radiation efficiency P_rad / P_acc within 5 % of 1 (lossless antenna)
   dipole: nf2ff P_rad (surface Poynting flux) within 2 % of the far-field integral
   normalised pattern |E(theta)| within 0.05 of the analytic pattern
   pattern independent of phi (rotational symmetry) within 1 %
   null on the antenna axis below -20 dB

 Tested with
  - python 3.13
  - openEMS v0.0.37+

 (c) 2026 Sean Mollet <sean@malmoset.com>

"""

import os, tempfile
import numpy as np

from CSXCAD  import ContinuousStructure
from openEMS import openEMS
from openEMS.physical_constants import C0, Z0
from openEMS.ports import CurvePort, LumpedPort

### Geometry
unit     = 1e-3   # drawing unit: mm
arm_len  = 75     # dipole arm / monopole length (mm)
mesh_res = 5.0    # mm

### Frequency
f_start = 0.7e9
f_stop  = 1.3e9
f_0     = 1.0e9
lambda0 = C0 / f_0 / unit    # 300 mm
pad     = lambda0/4 + 10*mesh_res   # clearance to the boundary, incl. PML

theta = np.arange(0, 181, 2.0)   # deg
phi   = np.array([0, 45, 90])     # deg


def analytic_pattern(theta_deg):
    th = np.deg2rad(theta_deg)
    on_axis = np.abs(np.sin(th)) < 1e-9
    with np.errstate(divide='ignore', invalid='ignore'):
        E = np.cos(np.pi/2*np.cos(th)) / np.sin(th)
    return np.where(on_axis, 0, E)


def setup(monopole):
    FDTD = openEMS(EndCriteria=1e-4)
    FDTD.SetGaussExcite(0.5*(f_start + f_stop), 0.5*(f_stop - f_start))
    if monopole:
        FDTD.SetBoundaryCond(['PML_8', 'PML_8', 'PML_8', 'PML_8', 'PEC', 'PML_8'])
    else:
        FDTD.SetBoundaryCond(['PML_8'] * 6)

    CSX = ContinuousStructure()
    FDTD.SetCSX(CSX)
    mesh = CSX.GetGrid()
    mesh.SetDeltaUnit(unit)
    mesh.AddLine('x', [-pad, 0, pad])
    mesh.AddLine('y', [-pad, 0, pad])
    if monopole:
        mesh.AddLine('z', [0, mesh_res, arm_len + pad])
    else:
        mesh.AddLine('z', [-(arm_len + pad), 0, arm_len + pad])
    mesh.SmoothMeshLines('all', mesh_res, ratio=1.4)
    return FDTD, CSX


def check_far_field(name, port, nf2ff, Sim_Path, D_ana, theta_range, nf2ff_theta, f_rec=None, half_space=False):
    """ Evaluate the far field at the S11 minimum, or at f_rec for a FD recording box.

    half_space: take P_rad and Dmax from the far-field integral over theta = 0..90 deg
    (the nf2ff P_rad / Dmax of a mirrored box include the image half space).
    """
    freq = np.linspace(f_start, f_stop, 201)
    port.CalcPort(Sim_Path, freq)
    s11_dB = 20*np.log10(np.abs(port.uf_ref / port.uf_inc))
    f_res  = freq[np.argmin(s11_dB)]
    print(f'  resonance at {f_res/1e9:.3f} GHz, S11 = {np.min(s11_dB):.1f} dB')
    if f_rec is not None:
        assert abs(f_res - f_rec) / f_rec < 0.05, \
            f'FAIL [{name}]: resonance at {f_res/1e9:.3f} GHz, expected near {f_rec/1e9:.3f} GHz (+/- 5 %)'
        f_res = f_rec

    port.CalcPort(Sim_Path, f_res)
    P_acc = port.P_acc[0]

    res = nf2ff.CalcNF2FF(Sim_Path, f_res, nf2ff_theta, phi)
    E   = res.E_norm[0]                   # (theta, phi) at r = 1 m
    E_n = E / np.max(E)

    # radiated power from the far-field radiation intensity U = |E|^2 r^2 / (2 Z0),
    # the pattern is rotationally symmetric: average over phi
    th   = np.deg2rad(nf2ff_theta)
    U    = E**2 / (2*Z0)
    y    = np.mean(U, axis=1) * np.sin(th)
    P_ff = 2*np.pi * np.sum(0.5*(y[1:] + y[:-1]) * np.diff(th))
    D_ff = 4*np.pi * np.max(U) / P_ff

    if half_space:
        Dmax, P_rad = D_ff, P_ff
    else:
        Dmax, P_rad = res.Dmax[0], res.Prad[0]
        print(f'  nf2ff P_rad / far-field integral = {res.Prad[0]/P_ff:.4f}')
        assert abs(res.Prad[0]/P_ff - 1) < 0.02, \
            f'FAIL [{name}]: nf2ff P_rad differs {abs(res.Prad[0]/P_ff - 1)*100:.1f} % from the far-field integral, expected < 2 %'

    eff = P_rad / P_acc
    print(f'  Dmax = {Dmax:.3f} ({10*np.log10(Dmax):.2f} dBi),  analytic {D_ana:.2f}')
    print(f'  P_rad / P_acc = {eff:.3f}')

    assert abs(Dmax - D_ana) / D_ana < 0.05, \
        f'FAIL [{name}]: Dmax = {Dmax:.3f}, expected {D_ana:.2f} (+/- 5 %)'
    assert abs(eff - 1) < 0.05, \
        f'FAIL [{name}]: P_rad / P_acc = {eff:.3f}, expected 1 (+/- 5 %)'

    mask = (nf2ff_theta >= theta_range[0]) & (nf2ff_theta <= theta_range[1])
    E_ana = analytic_pattern(nf2ff_theta)[:, None]
    pat_err = np.max(np.abs(E_n[mask] - E_ana[mask]))
    phi_err = np.max(np.abs(E_n[mask] - E_n[mask][:, :1]))
    axis_dB = 20*np.log10(E_n[0, 0] + 1e-12)
    print(f'  max pattern error = {pat_err:.3f},  phi asymmetry = {phi_err:.4f},  on-axis = {axis_dB:.1f} dB')

    assert pat_err < 0.05, \
        f'FAIL [{name}]: normalised pattern deviates {pat_err:.3f} from analytic, expected < 0.05'
    assert phi_err < 0.01, \
        f'FAIL [{name}]: pattern varies {phi_err:.4f} with phi, expected < 0.01 (rotational symmetry)'
    assert axis_dB < -20, \
        f'FAIL [{name}]: on-axis field {axis_dB:.1f} dB, expected a null < -20 dB'
    return res


### 1. Half-wave dipole, TD recording
print('Half-wave dipole (TD nf2ff recording)...')
Sim_Path = os.path.join(tempfile.gettempdir(), 'NF2FF_Dipole')
FDTD, CSX = setup(monopole=False)
port = CurvePort(CSX, 1, R=73, start=[0, 0, -arm_len], stop=[0, 0, arm_len], excite=1)
nf2ff = FDTD.CreateNF2FFBox()
FDTD.Run(Sim_Path, cleanup=True)
res_dipole = check_far_field('dipole', port, nf2ff, Sim_Path, 1.64, (0, 180), theta)
print('PASS [dipole]')

### 2. Monopole on a PEC ground plane, FD recording with PEC mirror
print('Quarter-wave monopole on ground (FD nf2ff recording, PEC mirror)...')
Sim_Path = os.path.join(tempfile.gettempdir(), 'NF2FF_Monopole')
FDTD, CSX = setup(monopole=True)
port = LumpedPort(CSX, 1, 36.5, [0, 0, 0], [0, 0, mesh_res], 'z', excite=1)
wire = CSX.AddMetal('monopole')
wire.AddCurve([[0, 0], [0, 0], [mesh_res, arm_len]], priority=10)
f_nf2ff = 0.91e9   # FD recording at the expected resonance
nf2ff = FDTD.CreateNF2FFBox(frequency=[f_nf2ff])
FDTD.Run(Sim_Path, cleanup=True)
res_mono = check_far_field('monopole', port, nf2ff, Sim_Path, 3.28, (0, 90), theta[theta <= 90],
                           f_rec=f_nf2ff, half_space=True)
print('PASS [monopole]')

print('PASS')

if 0:  # set to 1 for debugging plots
    import matplotlib.pyplot as plt

    fig = plt.figure(num='Far field', tight_layout=True)
    ax  = fig.add_subplot(111, polar=True)
    for res, lbl in ((res_dipole, 'dipole'), (res_mono, 'monopole')):
        E = res.E_norm[0][:, 0]
        ax.plot(res.theta, E/np.max(E), linewidth=2, label=lbl)   # res.theta is in rad
    ax.plot(np.deg2rad(theta), analytic_pattern(theta), 'k--', label='analytic')
    ax.set_theta_zero_location('N')
    ax.legend()
    plt.show()
