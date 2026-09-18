# -*- coding: utf-8 -*-
"""
 Steady-state detection test — sinusoidal excitation

 A sinusoidal excitation enables the steady-state detection, which ends the
 simulation once the field energy no longer changes from one excitation
 period to the next.  A TEM plane wave (Ex, travelling in +z) is launched in
 a 1D channel (PEC walls normal to x, PMC walls normal to y, PML at both
 z-ends) through a lossy dielectric slab.  The steady-state transmission,
 normalised to a reference run without the slab, is compared to the analytic
 slab transmission at the excitation frequency.

 Pass criteria
   both runs end by steady-state detection, long before the max. number
     of timesteps (end time < 50 excitation periods)
   the probe signal is periodic at the end: the last two periods differ < 0.1 %
   |T| within 0.1 dB and arg(T) within 2 deg of the analytic slab transmission

 Tested with
  - python 3.13
  - openEMS v0.0.37+

 (c) 2026 Sean Mollet <sean@malmoset.com>

"""

import os, tempfile
import numpy as np

from CSXCAD  import ContinuousStructure
from openEMS import openEMS
from openEMS.physical_constants import C0, EPS0

### Geometry
unit     = 1e-3    # drawing unit: mm
dz       = 0.5     # axial mesh resolution (mm)
width    = 1.0     # transverse channel width (mm)
length   = 150     # channel length (mm)
z_src    = 20      # excitation plane
slab_z0  = 60      # slab start
slab_d   = 20      # slab thickness (mm)
z_probe  = 110     # transmitted-field probe

### Excitation & slab material
f0      = 3e9
eps_r   = 4.0
kappa   = 0.05     # S/m
NrTS    = 1000000  # far more than needed, the steady-state detection must end the run
T0      = 1 / f0


def run(Sim_Path, slab=False):
    """ Run the channel and return the probe signal (t, Ex). """
    FDTD = openEMS(NrTS=NrTS, EndCriteria=1e-6)
    FDTD.SetSinusExcite(f0)
    FDTD.SetBoundaryCond(['PEC', 'PEC', 'PMC', 'PMC', 'PML_8', 'PML_8'])

    CSX = ContinuousStructure()
    FDTD.SetCSX(CSX)
    mesh = CSX.GetGrid()
    mesh.SetDeltaUnit(unit)
    mesh.AddLine('x', [0, width/2, width])
    mesh.AddLine('y', [0, width/2, width])
    mesh.AddLine('z', np.arange(0, length + dz, dz))

    exc = CSX.AddExcitation('plane', exc_type=0, exc_val=[1, 0, 0])
    exc.AddBox([0, 0, z_src], [width, width, z_src])

    CSX.AddProbe('et_trans', p_type=2).AddPoint([width/2, width/2, z_probe])

    if slab:
        mat = CSX.AddMaterial('slab', epsilon=eps_r, kappa=kappa)
        mat.AddBox([0, 0, slab_z0], [width, width, slab_z0 + slab_d], priority=10)

    FDTD.Run(Sim_Path, cleanup=True)

    # field probe file columns: t/s, Ex, Ey, Ez
    data = np.loadtxt(os.path.join(Sim_Path, 'et_trans'), comments='%')
    return data[:, 0], data[:, 1]


def phasor(t, val, t_start, t_stop):
    """ least-squares fit of a cos(w t) + b sin(w t) in [t_start, t_stop), returns a - jb """
    m = (t >= t_start) & (t < t_stop)
    w = 2*np.pi*f0
    A = np.column_stack((np.cos(w*t[m]), np.sin(w*t[m])))
    (a, b), *_ = np.linalg.lstsq(A, val[m], rcond=None)
    return a - 1j*b


def steady_state(label, t, val):
    t_end = t[-1]
    print(f'  {label}: simulation ended at {t_end*1e9:.2f} ns = {t_end/T0:.0f} periods')
    assert t_end < 50*T0, \
        f'FAIL [{label}]: simulation ran for {t_end/T0:.0f} periods, expected the steady-state detection to end it within 50'

    last = phasor(t, val, t_end - T0,   t_end)
    prev = phasor(t, val, t_end - 2*T0, t_end - T0)
    change = abs(last - prev) / abs(last)
    print(f'  {label}: change between the last two periods: {change*100:.4f} %')
    assert change < 1e-3, \
        f'FAIL [{label}]: probe changes {change*100:.3f} % between the last two periods, expected < 0.1 % (not steady)'
    return phasor(t, val, t_end - 4*T0, t_end)


print('Running reference (no slab)...')
t_ref, E_ref_t = run(os.path.join(tempfile.gettempdir(), 'SteadyState_Ref'))
E_ref = steady_state('reference', t_ref, E_ref_t)

print('Running lossy dielectric slab...')
t_slab, E_slab_t = run(os.path.join(tempfile.gettempdir(), 'SteadyState_Slab'), slab=True)
E_slab = steady_state('slab', t_slab, E_slab_t)

### analytic slab transmission at f0 (e^{jwt} convention)
w      = 2*np.pi*f0
eps_c  = eps_r - 1j*kappa/(w*EPS0)
k0     = w / C0
k      = k0*np.sqrt(eps_c)
k      = -k if np.imag(k) > 0 else k
eta    = 1/np.sqrt(eps_c)
d      = slab_d*unit
T_ana  = np.exp(1j*k0*d) / (np.cos(k*d) + 0.5j*(eta + 1/eta)*np.sin(k*d))

T_fdtd = E_slab / E_ref
dB_err    = 20*np.log10(abs(T_fdtd/T_ana))
phase_err = np.angle(T_fdtd/T_ana, deg=True)
print(f'T = {20*np.log10(abs(T_fdtd)):.2f} dB, {np.angle(T_fdtd, deg=True):.1f} deg  '
      f'(analytic {20*np.log10(abs(T_ana)):.2f} dB, {np.angle(T_ana, deg=True):.1f} deg)')

assert abs(dB_err) < 0.1, \
    f'FAIL: |T| deviates {dB_err:.3f} dB from analytic, expected < 0.1 dB'
assert abs(phase_err) < 2, \
    f'FAIL: arg(T) deviates {phase_err:.2f} deg from analytic, expected < 2 deg'

print('PASS')

if 0:  # set to 1 for debugging plots
    import matplotlib.pyplot as plt

    fig, axis = plt.subplots(num='Steady state', tight_layout=True)
    axis.plot(t_ref*1e9,  E_ref_t,  'k-',  linewidth=1, label='reference')
    axis.plot(t_slab*1e9, E_slab_t, 'r-',  linewidth=1, label='slab')
    axis.grid()
    axis.set_xlabel('Time (ns)')
    axis.set_ylabel('Ex (V/m)')
    axis.legend()
    plt.show()
