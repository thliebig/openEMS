# -*- coding: utf-8 -*-
"""
 Multi-pole Debye test -- FR4-like substrate, propagation constant and slab

 The companion to Debye_Water.py.  Where that one stresses the Debye ADE with a
 strong pole (d_eps / eps_inf = 14.2), this one covers the case the ADE has
 actually been used for since it was added: a low-loss substrate with a flat
 loss tangent over a wide band, modelled as a sum of weak Debye poles
 (Djordjevic-Sarkar).  Five poles with equal d_eps, relaxation frequencies
 log-spaced 0.1 .. 50 GHz:

     eps_r(f) = eps_inf + sum_k d_eps / (1 + j w tau_k)
     eps_inf = 4.10,  d_eps = 0.09 per pole,  f_relax = 0.1, 0.5, 2.5, 12.5, 50 GHz

 which realises eps' = 4.25 .. 4.38 and tand = 0.0194 .. 0.0207 over 1 .. 10 GHz
 -- FR4-like, and flat to +-3% across a decade.  sum(d_eps)/eps_inf = 0.11, so
 this sits far inside the stable range of the ADE and must pass before and after
 any change to it: it is the regression guard for the working case, and the only
 coverage of the multi-pole (order > 1) path.

 Two measurements, in that order of importance.

 1. Propagation constant.  Two probes inside a long block of the material,
    50 mm apart, give the ratio

        E2/E1 = exp(-alpha d) exp(-j beta d)

    so beta and alpha come out directly, with no interface and no reference run
    in the way.  The record is gated before the block's far face echoes back to
    the second probe -- the analytic phase velocity is used to schedule that
    window, not to produce the result.  The probe spacing is 50 mm so that
    beta*d < pi at 1 GHz and the frequency unwrap is unambiguous without being
    anchored to the analytic answer.  Compared against

        k = k0 sqrt(eps_r),  beta = Re k,  alpha = -Im k    (e^{jwt}, Im k <= 0)

 2. Slab transmission, which the first measurement deliberately excludes: a
    30 mm slab, normalised to an empty reference run, against

        T(f) = exp(j k0 d) / (cos(k d) + j/2 (eta + 1/eta) sin(k d))

    with eta = 1/sqrt(eps_r).  30 mm puts four Fabry-Perot resonances in band
    (2.4, 4.8, 7.2, 9.7 GHz), so this is sensitive to eps' through the resonance
    frequencies and to tand through their depth, and it covers the staircased
    eps_r / vacuum interfaces.

 The alpha criterion is the sharp one.  Comparing the z-domain transfer function
 of the ADE against the analytic permittivity suggests this material is
 insensitive to how the Debye branch is discretised -- the dispersive part is
 only 11% of eps_r, and the two agree to better than 0.001%.  Measurement says
 otherwise: alpha moves by a factor of 3 between two discretisations of the same
 branch, so that comparison misses something, most likely the half-timestep
 offset with which the branch current enters the voltage update.  Trust alpha
 here, not the transfer function.

 Pass criteria (1 .. 10 GHz)
   max relative error of beta  < 0.02%
   max relative error of alpha < 1%
   max |T_fdtd - T_analytic|   < 0.05

 This test was written together with the Debye ADE rework, not before it, and
 the tolerances are set to what that rework delivers with a factor of three in
 hand: beta 0.006%, alpha 0.32%, |dT| 0.0045 at dz = 0.1 mm.  The scheme it
 replaced -- which integrated the branch's algebraic loop equation rather than
 the capacitor equation -- gave beta 0.027% and alpha 2.05% on the same mesh, so
 alpha is the criterion that would catch a regression to it.

 alpha converges slowly with the mesh because it is a small quantity; the
 underlying amplitude error over the 50 mm baseline is a few parts in a thousand.
 Do not coarsen dz without re-measuring it.

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

### Mesh
unit     = 1e-3    # drawing unit: mm
dz       = 0.1     # axial mesh resolution (mm)
width    = 1.0     # transverse channel width (mm)

### Frequency
f_start = 1e9
f_stop  = 10e9
freq    = np.logspace(np.log10(f_start), np.log10(f_stop), 141)
w       = 2*np.pi*freq

### FR4-like 5-pole Djordjevic-Sarkar set
eps_inf = 4.10
d_eps   = 0.09
f_relax = np.array([0.1, 0.5, 2.5, 12.5, 50.0]) * 1e9
tau     = 1 / (2*np.pi*f_relax)

### The timestep openEMS will pick: uniform cartesian mesh, vacuum limited.
### It matters here: openEMS silently drops any pole with tau <= 2 dT -- see the
### "relaxation time is to small" guard in
### Operator_Ext_LorentzMaterial::BuildExtension() -- so a coarser mesh would
### quietly turn this into a four-pole material.
dx = dy = width/2 * unit
dT = 1 / (C0 * np.sqrt(1/dx**2 + 1/dy**2 + 1/(dz*unit)**2))


def eps_substrate(f):
    return eps_inf + sum(d_eps / (1 + 1j*2*np.pi*f*t) for t in tau)


def add_material(CSX, z0, z1):
    mat = CSPropDebyeMaterial(CSX.GetParameterSet(), order=len(tau), epsilon=eps_inf)
    mat.SetName('fr4')
    for k, t in enumerate(tau):
        mat.SetDispersiveMaterialProperty(k, eps_delta=d_eps, eps_relax=t)
    CSX.AddProperty(mat)
    mat.AddBox([0, 0, z0], [width, width, z1], priority=10)
    return mat


def run(Sim_Path, length, z_src, probes, material=None):
    """ 1D channel of `length` mm, optional material block (z0, z1), probes at
        the given z positions.  Returns t and one Ex column per probe. """
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

    for n, z in enumerate(probes):
        CSX.AddProbe('et_{}'.format(n), p_type=2).AddPoint([width/2, width/2, z])

    if material is not None:
        add_material(CSX, *material)

    FDTD.Run(Sim_Path, cleanup=True, exact_endcriteria=True)

    # field probe file columns: t/s, Ex, Ey, Ez
    out = [np.loadtxt(os.path.join(Sim_Path, 'et_{}'.format(n)), comments='%')
           for n in range(len(probes))]
    return out[0][:, 0], [d[:, 1] for d in out]


def slab_transmission(eps_r, d_mm):
    k0  = w / C0
    k   = k0 * np.sqrt(eps_r + 0j)
    k   = np.where(np.imag(k) > 0, -k, k)   # passive medium: Im(k) <= 0 for e^{jwt}
    eta = np.sqrt(1 / eps_r + 0j)
    eta = np.where(np.real(eta) < 0, -eta, eta)
    d   = d_mm * unit
    return np.exp(1j*k0*d) / (np.cos(k*d) + 0.5j*(eta + 1/eta)*np.sin(k*d))


eps_r = eps_substrate(freq)
tand  = -np.imag(eps_r) / np.real(eps_r)
k_ana = w/C0 * np.sqrt(eps_r + 0j)
k_ana = np.where(np.imag(k_ana) > 0, -k_ana, k_ana)
beta_ana, alpha_ana = np.real(k_ana), -np.imag(k_ana)

print('{} poles, sum(d_eps)/eps_inf = {:.3f}'.format(len(tau), len(tau)*d_eps/eps_inf))
print("eps' = {:.3f} .. {:.3f},  tand = {:.5f} .. {:.5f}  over {:.0f} .. {:.0f} GHz".format(
      eps_r.real.min(), eps_r.real.max(), tand.min(), tand.max(),
      f_start/1e9, f_stop/1e9))
print('fastest pole: tau = {:.2f} ps, dT = {:.3f} ps, tau/dT = {:.1f}'.format(
      tau.min()*1e12, dT*1e12, tau.min()/dT))
assert tau.min() > 2*dT, \
    'FAIL: fastest pole tau = {:.2f} ps <= 2 dT = {:.2f} ps -- openEMS would drop ' \
    'it and simulate a {}-pole material instead'.format(
        tau.min()*1e12, 2*dT*1e12, len(tau)-1)


### 1. propagation constant, from two probes inside a long block
blk_z0, blk_z1 = 30, 250
z_src          = 10
z_p1, z_p2     = 80, 130
d_probe        = (z_p2 - z_p1) * unit

print('\nRunning block (alpha, beta)...')
t, (E1_t, E2_t) = run(os.path.join(tempfile.gettempdir(), 'Debye_Substrate_Block'),
                      280, z_src, [z_p1, z_p2], material=(blk_z0, blk_z1))

### Gate before the far face of the block echoes back to the second probe. The
### analytic phase velocity only schedules the window; the result comes from the
### probes.
v_mid  = C0 / np.sqrt(eps_substrate(np.sqrt(f_start*f_stop)).real)
t_echo = ((blk_z0 - z_src)*unit/C0
          + (z_p2 - blk_z0)*unit/v_mid
          + 2*(blk_z1 - z_p2)*unit/v_mid)
gate = t <= t_echo
print('record {} samples, gate keeps {} (echo at p2 by {:.2f} ns)'.format(
      len(t), int(gate.sum()), t_echo*1e9))
assert gate.sum() > 50, 'FAIL: gate kept only {} samples'.format(int(gate.sum()))

E1 = DFT_time2freq(t[gate], E1_t[gate], freq)
E2 = DFT_time2freq(t[gate], E2_t[gate], freq)
ratio      = E2 / E1
alpha_fdtd = -np.log(np.abs(ratio)) / d_probe
beta_fdtd  = np.unwrap(-np.angle(ratio)) / d_probe

assert beta_ana[0]*d_probe < np.pi, \
    'FAIL: beta*d = {:.2f} rad at {:.1f} GHz, the phase unwrap is ambiguous -- ' \
    'move the probes closer together'.format(beta_ana[0]*d_probe, f_start/1e9)

err_beta  = np.max(np.abs(beta_fdtd/beta_ana - 1))
err_alpha = np.max(np.abs(alpha_fdtd/alpha_ana - 1))
print('beta : {:.1f} .. {:.1f} rad/m, max rel error {:.4f}%'.format(
      beta_fdtd.min(), beta_fdtd.max(), err_beta*100))
print('alpha: {:.3f} .. {:.3f} Np/m, max rel error {:.2f}%'.format(
      alpha_fdtd.min(), alpha_fdtd.max(), err_alpha*100))
assert err_beta < 2e-4, \
    'FAIL: max relative error of beta = {:.4f}%, expected < 0.02%'.format(err_beta*100)
assert err_alpha < 1e-2, \
    'FAIL: max relative error of alpha = {:.2f}%, expected < 1%'.format(err_alpha*100)


### 2. slab transmission, which also covers the interfaces
slab_z0, slab_d, z_probe = 30, 30.0, 80

print('\nRunning reference (no slab)...')
t_ref, (E_ref_t,) = run(os.path.join(tempfile.gettempdir(), 'Debye_Substrate_Ref'),
                        100, z_src, [z_probe])
print('Running slab...')
t_slab, (E_slab_t,) = run(os.path.join(tempfile.gettempdir(), 'Debye_Substrate_Slab'),
                          100, z_src, [z_probe],
                          material=(slab_z0, slab_z0 + slab_d))

ref_peak = np.max(np.abs(E_ref_t))
finite   = bool(np.all(np.isfinite(E_slab_t)))
peak     = np.max(np.abs(E_slab_t)) if finite else np.inf
print('max |Ex| behind the slab: {:.4g}  (reference {:.4g})'.format(peak, ref_peak))
assert finite and peak < 10 * ref_peak, \
    'FAIL: the Debye slab diverged (max |Ex| = {:.4g})'.format(peak)
assert peak > 0.01 * ref_peak, \
    'FAIL: no transmitted pulse behind the Debye slab (max |Ex| = {:.4g})'.format(peak)

T_fdtd = DFT_time2freq(t_slab, E_slab_t, freq) / DFT_time2freq(t_ref, E_ref_t, freq)
T_ana  = slab_transmission(eps_r, slab_d)
err_T  = np.max(np.abs(T_fdtd - T_ana))
print('max |T_fdtd - T_analytic| = {:.4f}'.format(err_T))
assert err_T < 0.05, \
    'FAIL: max |T_fdtd - T_analytic| = {:.4f}, expected < 0.05'.format(err_T)

print('PASS')

if 0:  # set to 1 for debugging plots
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(3, 1, num='FR4-like substrate',
                             tight_layout=True, sharex=True)
    axes[0].semilogx(freq/1e9, beta_fdtd, 'C0-',  linewidth=2, label='openEMS')
    axes[0].semilogx(freq/1e9, beta_ana,  'C0--', linewidth=1, label='analytic')
    axes[1].semilogx(freq/1e9, alpha_fdtd, 'C1-',  linewidth=2)
    axes[1].semilogx(freq/1e9, alpha_ana,  'C1--', linewidth=1)
    axes[2].semilogx(freq/1e9, 20*np.log10(np.abs(T_fdtd)), 'C2-',  linewidth=2)
    axes[2].semilogx(freq/1e9, 20*np.log10(np.abs(T_ana)),  'C2--', linewidth=1)
    axes[0].set_ylabel('beta (rad/m)')
    axes[0].legend()
    axes[1].set_ylabel('alpha (Np/m)')
    axes[2].set_ylabel('|T| 30 mm slab (dB)')
    axes[2].set_xlabel('Frequency (GHz)')
    for ax in axes:
        ax.grid(which='both')
    plt.show()
