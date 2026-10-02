# -*- coding: utf-8 -*-
"""
 Debye material test -- plane-wave transmission through a water slab

 The same 1D channel as Dispersive_Materials.py, but with a Debye material of
 realistic strength: water at 25 C as a single Debye pole,

     eps_r(f) = eps_inf + d_eps / (1 + j w tau)
     eps_inf = 5.16,  d_eps = 73.2  (eps_static = 78.36),  tau = 8.27 ps

 so d_eps / eps_inf = 14.2.  This is the case that used to diverge: before the
 Debye ADE was reworked it was unstable above a ratio of roughly 1, which excluded
 most physical Debye media, and it failed silently -- the field went to NaN, the
 NaN energy satisfied the end criteria, and the run stopped reporting success.
 The test was written together with that rework, not before it.

 The transmitted field, normalised to a reference run without the slab, is
 compared to the analytic slab transmission

     T(f) = exp(j k0 d) / (cos(k d) + j/2 (eta + 1/eta) sin(k d))

 with k = k0 sqrt(eps_r) and eta = 1/sqrt(eps_r)  (e^{jwt} convention).

 Pass criteria (2 .. 9 GHz)
   the fields stay finite
   max |T_fdtd - T_analytic| < 0.05

 Mesh: lambda inside the water is ~4 mm at 9 GHz (|sqrt(eps_r)| ~ 8.4), so
 dz = 0.05 mm is lambda/80 and the 2 mm slab is resolved by 40 cells.  That is
 not generous for its own sake: the error is dominated by the staircasing of the
 eps_r ~ 78 / 1 interfaces and falls about linearly with dz, measured

     dz / mm      0.2     0.1     0.05
     max |dT|     --      0.053   0.026

 so dz = 0.1 mm does not clear 0.05 and leaves no margin at all.

 Tested with
  - python 3.13
  - openEMS v0.0.37+

"""

import os, tempfile
import numpy as np

from CSXCAD  import ContinuousStructure
from CSXCAD.CSProperties import CSPropDebyeMaterial
from openEMS import openEMS
from openEMS.physical_constants import C0
from openEMS.utilities import DFT_time2freq

### Geometry
unit     = 1e-3    # drawing unit: mm
dz       = 0.05    # axial mesh resolution (mm)
width    = 1.0     # transverse channel width (mm)
length   = 80      # channel length (mm)
z_src    = 10      # excitation plane
slab_z0  = 30      # slab start
slab_d   = 2.0     # slab thickness (mm)
z_probe  = 50      # transmitted-field probe

### Frequency
f_start = 1e9
f_stop  = 10e9
freq    = np.linspace(2e9, 9e9, 141)   # evaluation band
w       = 2*np.pi*freq

### Water at 25 C, single Debye pole
eps_inf = 5.16
d_eps   = 73.2     # eps_static = 78.36
tau     = 8.27e-12

### The timestep openEMS will pick: uniform cartesian mesh, vacuum limited.
### Needed here because the Debye pole is silently dropped, and the material
### then behaves as plain eps_inf, when tau <= 2 dT -- see the "relaxation time
### is to small" guard in Operator_Ext_LorentzMaterial::BuildExtension().
dx = dy = width/2 * unit
dT = 1 / (C0 * np.sqrt(1/dx**2 + 1/dy**2 + 1/(dz*unit)**2))


def run(Sim_Path, with_slab=False):
    """ Run the channel with or without the water slab, return Ex(freq). """
    FDTD = openEMS(NrTS=100000, EndCriteria=1e-6)
    FDTD.SetGaussExcite(0.5*(f_start + f_stop), 0.5*(f_stop - f_start))
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

    probe = CSX.AddProbe('et_trans', p_type=2)
    probe.AddPoint([width/2, width/2, z_probe])

    if with_slab:
        mat = CSPropDebyeMaterial(CSX.GetParameterSet(), order=1, epsilon=eps_inf)
        mat.SetName('water')
        mat.SetDispersiveMaterialProperty(0, eps_delta=d_eps, eps_relax=tau)
        CSX.AddProperty(mat)
        mat.AddBox([0, 0, slab_z0], [width, width, slab_z0 + slab_d], priority=10)

    FDTD.Run(Sim_Path, cleanup=True, exact_endcriteria=True)

    # field probe file columns: t/s, Ex, Ey, Ez
    data = np.loadtxt(os.path.join(Sim_Path, 'et_trans'), comments='%')
    return data[:, 0], data[:, 1]


def slab_transmission(eps_r):
    k0  = w / C0
    k   = k0 * np.sqrt(eps_r + 0j)
    k   = np.where(np.imag(k) > 0, -k, k)   # passive medium: Im(k) <= 0 for e^{jwt}
    eta = np.sqrt(1 / eps_r + 0j)
    eta = np.where(np.real(eta) < 0, -eta, eta)
    d   = slab_d * unit
    return np.exp(1j*k0*d) / (np.cos(k*d) + 0.5j*(eta + 1/eta)*np.sin(k*d))


def eps_water(f):
    return eps_inf + d_eps / (1 + 1j*2*np.pi*f*tau)


print('d_eps / eps_inf = {:.1f}'.format(d_eps/eps_inf))
print('tau = {:.2f} ps, dT = {:.3f} ps, tau/dT = {:.1f}'.format(
      tau*1e12, dT*1e12, tau/dT))
assert tau > 2*dT, \
    'FAIL: tau = {:.2f} ps <= 2 dT = {:.2f} ps -- openEMS would drop the Debye ' \
    'pole and simulate plain eps_inf instead'.format(tau*1e12, 2*dT*1e12)

print('Running reference (no slab)...')
t_ref, E_ref_t = run(os.path.join(tempfile.gettempdir(), 'Debye_Water_Ref'))

print('Running water slab...')
t_slab, E_slab_t = run(os.path.join(tempfile.gettempdir(), 'Debye_Water'), True)

### Divergence guards.  A diverging Debye slab does not announce itself: the
### field goes to NaN, the NaN energy compares as below the end criteria, and
### openEMS reports "end-criteria ... reached (-nandB)" and stops -- a run that
### looks finished, with a probe record of zeros because it quit before the
### pulse arrived.  So check the transmitted pulse is there and is bounded,
### not just that the numbers are finite.
ref_peak = np.max(np.abs(E_ref_t))
finite   = bool(np.all(np.isfinite(E_slab_t)))
peak     = np.max(np.abs(E_slab_t)) if finite else np.inf
print('record length: {} samples (reference {})'.format(len(E_slab_t), len(E_ref_t)))
print('max |Ex| behind the slab: {:.4g}  (reference {:.4g})'.format(peak, ref_peak))

assert finite, \
    'FAIL: the field behind the Debye slab is not finite -- the slab diverged'
assert peak < 10 * ref_peak, \
    'FAIL: the Debye slab diverged, max |Ex| = {:.4g} is {:.3g}x the incident ' \
    'field'.format(peak, peak/ref_peak)
assert peak > 0.01 * ref_peak, \
    'FAIL: no transmitted pulse behind the Debye slab (max |Ex| = {:.4g}, ' \
    '{:.3g}x the incident field). If the record is short, the run diverged and ' \
    'stopped on a NaN end criteria -- look for "-nandB" in the log.'.format(
        peak, peak/ref_peak)

T_fdtd = DFT_time2freq(t_slab, E_slab_t, freq) / DFT_time2freq(t_ref, E_ref_t, freq)
T_ana  = slab_transmission(eps_water(freq))

err = np.max(np.abs(T_fdtd - T_ana))
print('max |T_fdtd - T_analytic| = {:.4f}'.format(err))
assert err < 0.05, \
    'FAIL: max |T_fdtd - T_analytic| = {:.4f}, expected < 0.05'.format(err)

print('PASS')

if 0:  # set to 1 for debugging plots
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 1, num='Water slab transmission',
                             tight_layout=True, sharex=True)
    axes[0].plot(freq/1e9, 20*np.log10(np.abs(T_fdtd)), 'C0-', linewidth=2,
                 label='openEMS')
    axes[0].plot(freq/1e9, 20*np.log10(np.abs(T_ana)), 'C0--', linewidth=1,
                 label='analytic')
    axes[1].plot(freq/1e9, np.angle(T_fdtd, deg=True), 'C0-', linewidth=2)
    axes[1].plot(freq/1e9, np.angle(T_ana,  deg=True), 'C0--', linewidth=1)
    axes[0].set_ylabel('|T| (dB)')
    axes[0].legend()
    axes[1].set_ylabel('arg(T) (deg)')
    axes[1].set_xlabel('Frequency (GHz)')
    for ax in axes:
        ax.grid()
    plt.show()
