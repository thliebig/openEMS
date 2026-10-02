# -*- coding: utf-8 -*-
"""
Modelling a flat-loss substrate with Debye poles, and what it costs per decade.

A real laminate like FR4 is specified as "eps_r = 4.3, tanD = 0.02" with both
roughly constant over the datasheet range.  openEMS has no constant-tanD
material, so the question is what to put in the model instead.  The obvious
choice, a plain material with a conductivity, is wrong in a way that is easy to
miss: a fixed `kappa` gives

    tanD(f) = kappa / (2 pi f eps0 eps_r)

which falls as 1/f.  Set it to 0.02 at 3 GHz and it is 0.04 at 1.5 GHz and
0.013 at 6 GHz -- a factor of three across one decade.  The attenuation of a
line built from it is then almost frequency independent, where a real flat-tanD
laminate loses roughly proportionally to frequency.

The fix is a sum of Debye poles.  This example fits three of them -- one, two
and five poles -- to a flat tanD, runs a 50 Ohm microstrip on each, and extracts
eps_eff and the attenuation back out of the line, so you can see directly how
far each fit holds.

  narrow   1 pole,  relaxation at  3 GHz
  medium   2 poles, relaxation 1 .. 10 GHz
  wide     5 poles, relaxation 0.1 .. 50 GHz

The reference is the Djordjevic-Sarkar model, the causal realisation of a
constant tanD over a finite band; a strictly constant tanD at all frequencies is
not causal and cannot be realised by any material.

What to expect: panels 1 and 2 are analytic and span 1 MHz .. 100 GHz, so the
bandwidth of each fit is visible.  Panels 3 and 4 are the simulated line over
1 .. 5 GHz, against the closed-form microstrip prediction

    eps_eff = (eps_r+1)/2 + (eps_r-1)/2 (1 + 12 h/w)^-1/2          (Hammerstad)
    alpha_d = 27.3 (eps_r/sqrt(eps_eff)) (eps_eff-1)/(eps_r-1) tanD/lambda0

Those closed forms are good to a few percent at best -- Hammerstad is static and
carries none of the microstrip dispersion the simulation shows -- so treat the
agreement in panels 3 and 4 as a sanity check, not as a measurement of either.
Measured here: alpha within about 2% of the formula at mid band, drifting to
about 10% at 4.5 GHz, and eps_eff 1.6% above Hammerstad at 1.5 GHz rising to 4%
at 4.5 GHz -- that rise is the dispersion the static formula does not have.

Three simulations of about 434k cells, roughly 80 s in total.

Tested with
 - python 3.13
 - openEMS v0.0.37+
"""

import os, tempfile
import numpy as np
import matplotlib.pyplot as plt

from CSXCAD import ContinuousStructure
from CSXCAD.CSProperties import CSPropDebyeMaterial
from openEMS import openEMS
from openEMS.physical_constants import C0, EPS0

post_proc_only = False      # True: skip the simulations and re-plot

### Target laminate -----------------------------------------------------------
EPS_TARGET  = 4.30          # eps_r, flat
TAND_TARGET = 0.020         # tanD, flat
DS_BAND     = (1e3, 1e12)   # Djordjevic-Sarkar support, Hz

### Microstrip ----------------------------------------------------------------
unit     = 1e-3
msl_w    = 2.9              # mm, ~50 Ohm on 1.5 mm FR4
sub_h    = 1.5              # mm
sub_w    = 40.0             # mm
msl_len  = 100.0            # mm
f0, fc   = 3e9, 2.5e9       # Gaussian excitation
freq     = np.linspace(1e9, 5e9, 201)
freq_ana = np.logspace(6, 11, 400)


def eps_djordjevic_sarkar(f, eps_r=EPS_TARGET, tand=TAND_TARGET,
                          band=DS_BAND, f_ref=3e9):
    """Causal constant-tanD reference,

        eps(w) = eps_inf + d_eps ln((1 + jw/w1)/(1 + jw/w2)) / ln(w2/w1)

    with d_eps and eps_inf set so that eps' = eps_r and tanD = tand at f_ref."""
    w1, w2 = 2*np.pi*band[0], 2*np.pi*band[1]
    kern = lambda w: np.log((1 + 1j*w/w1)/(1 + 1j*w/w2))/np.log(w2/w1)
    k_ref = kern(2*np.pi*f_ref)
    d_eps = eps_r*tand/(-k_ref.imag)
    eps_inf = eps_r - d_eps*k_ref.real
    return eps_inf + d_eps*kern(2*np.pi*f)


def fit_poles(n_poles, f_lo, f_hi, f_chk):
    """Equal-d_eps Debye poles, relaxation frequencies log-spaced over
       [f_lo, f_hi], with eps_inf and the common d_eps least-squared against the
       flat target over f_chk.  Returns (eps_inf, d_eps, tau[])."""
    tau = 1/(2*np.pi*np.logspace(np.log10(f_lo), np.log10(f_hi), n_poles))
    w = 2*np.pi*f_chk
    A_re = np.column_stack([np.ones_like(w), sum(1/(1+(w*t)**2)    for t in tau)])
    A_im = np.column_stack([np.zeros_like(w), sum(w*t/(1+(w*t)**2) for t in tau)])
    A = np.vstack([A_re, A_im])
    b = np.concatenate([np.full_like(w, EPS_TARGET),
                        np.full_like(w, EPS_TARGET*TAND_TARGET)])
    eps_inf, d_eps = np.linalg.lstsq(A, b, rcond=None)[0]
    return eps_inf, d_eps, tau


def eps_poles(f, fit):
    eps_inf, d_eps, tau = fit
    return eps_inf + sum(d_eps/(1 + 1j*2*np.pi*f*t) for t in tau)


### The three fits ------------------------------------------------------------
f_band = np.logspace(9, np.log10(5e9), 60)          # the simulated band
FITS = {
    'narrow (1 pole)': fit_poles(1,   3e9,  3e9, f_band),
    'medium (2 pole)': fit_poles(2,   1e9, 10e9, f_band),
    'wide (5 pole)'  : fit_poles(5, 0.1e9, 50e9, f_band),
}


def hammerstad(eps_r, w=msl_w, h=sub_h):
    return (eps_r+1)/2 + (eps_r-1)/2*(1 + 12*h/w)**-0.5


def alpha_d_dB_m(eps_r, eps_eff, tand, f):
    return 27.3*(eps_r/np.sqrt(eps_eff))*((eps_eff-1)/(eps_r-1))*tand/(C0/f)


def run_msl(tag, fit):
    """50 Ohm microstrip on the given pole set; returns the complex beta."""
    Sim_Path = os.path.join(tempfile.gettempdir(), 'MSL_Debye_' + tag)
    if True:   # the port object is needed for CalcPort(), so always build it
        F = openEMS(NrTS=100000, EndCriteria=1e-5)
        F.SetGaussExcite(f0, fc)
        # PEC under the ground plane, PML elsewhere
        F.SetBoundaryCond(['PML_8','PML_8','PML_8','PML_8','PEC','PML_8'])
        CSX = ContinuousStructure(); F.SetCSX(CSX)
        mesh = CSX.GetGrid(); mesh.SetDeltaUnit(unit)
        mesh.AddLine('x', [-sub_w/2, -msl_w/2, 0, msl_w/2, sub_w/2])
        mesh.AddLine('y', [0, msl_len])
        mesh.AddLine('z', [0, sub_h, sub_h+25])
        mesh.SmoothMeshLines('x', 0.5); mesh.SmoothMeshLines('y', 0.5)
        mesh.AddLine('z', np.linspace(0, sub_h, 7)); mesh.SmoothMeshLines('z', 1.5)

        eps_inf, d_eps, tau = fit
        sub = CSPropDebyeMaterial(CSX.GetParameterSet(), order=len(tau),
                                  epsilon=eps_inf)
        sub.SetName('substrate'); CSX.AddProperty(sub)
        for k, t in enumerate(tau):
            sub.SetDispersiveMaterialProperty(k, eps_delta=d_eps,
                                                 eps_relax=float(t))
        sub.AddBox([-sub_w/2, 0, 0], [sub_w/2, msl_len, sub_h], priority=2)

        pec = CSX.AddMetal('pec')
        pec.AddBox([-sub_w/2, 0, 0], [sub_w/2, msl_len, 0], priority=10)  # ground
        port = F.AddMSLPort(1, pec, [-msl_w/2, 0, sub_h], [msl_w/2, msl_len, 0],
                            'y', 'z', excite=-1, FeedShift=10,
                            MeasPlaneShift=msl_len/3, priority=10)
        if not post_proc_only:
            F.Run(Sim_Path, cleanup=True, verbose=0)
    port.CalcPort(Sim_Path, freq, ref_impedance=50)
    return port.beta


### Simulate ------------------------------------------------------------------
results = {}
for name, fit in FITS.items():
    print('Running {} ...'.format(name))
    beta = run_msl(name.split()[0], fit)
    results[name] = dict(eps_eff=(beta.real*C0/(2*np.pi*freq))**2,
                         alpha_dB=-beta.imag*8.686)

### Report --------------------------------------------------------------------
print('\n{:17s} {:>9} {:>9} {:>9} {:>9}'.format(
      '', 'eps_inf', 'd_eps_k', 'poles', 'sum/eps_inf'))
for name, (eps_inf, d_eps, tau) in FITS.items():
    print('{:17s} {:9.4f} {:9.5f} {:9d} {:9.4f}'.format(
          name, eps_inf, d_eps, len(tau), len(tau)*d_eps/eps_inf))

print('\ntanD over 1 .. 5 GHz (target {:.3f}):'.format(TAND_TARGET))
for name, fit in FITS.items():
    e = eps_poles(f_band, fit); td = -e.imag/e.real
    print('  {:17s} {:.4f} .. {:.4f}   ripple {:+.1f}%'.format(
          name, td.min(), td.max(), (td.max()/td.min()-1)*100))

print('\nmicrostrip, simulated vs closed form:')
for name, fit in FITS.items():
    r = results[name]
    for f in (1.5e9, 3e9, 4.5e9):
        i = np.argmin(abs(freq-f))
        e = eps_poles(f, fit)
        a_ana = alpha_d_dB_m(e.real, r['eps_eff'][i], -e.imag/e.real, f)
        print('  {:17s} {:.1f} GHz  eps_eff {:.4f} (Hammerstad {:.4f})  '
              'alpha {:6.2f} dB/m (formula {:6.2f}, {:+.1f}%)'.format(
              name, f/1e9, r['eps_eff'][i], hammerstad(e.real),
              r['alpha_dB'][i], a_ana, (r['alpha_dB'][i]/a_ana-1)*100))

### Plot ----------------------------------------------------------------------
fig, ax = plt.subplots(2, 2, figsize=(11, 7.5), tight_layout=True)
ds = eps_djordjevic_sarkar(freq_ana)

ax[0,0].semilogx(freq_ana/1e9, ds.real, 'k--', lw=1.5, label='Djordjevic-Sarkar')
ax[0,1].semilogx(freq_ana/1e9, -ds.imag/ds.real, 'k--', lw=1.5, label='Djordjevic-Sarkar')
for n, (name, fit) in enumerate(FITS.items()):
    e = eps_poles(freq_ana, fit)
    ax[0,0].semilogx(freq_ana/1e9, e.real, f'C{n}-', label=name)
    ax[0,1].semilogx(freq_ana/1e9, -e.imag/e.real, f'C{n}-', label=name)
for a in (ax[0,0], ax[0,1]):
    a.axvspan(1, 5, color='0.9', zorder=0)
    a.set_xlabel('f (GHz)'); a.grid(which='both'); a.legend(fontsize=8)
ax[0,0].set_ylabel("eps_r'"); ax[0,0].set_ylim(3.8, 5.2)
ax[0,1].set_ylabel('tan(delta)'); ax[0,1].set_ylim(0, 0.05)
ax[0,0].set_title('material model (shaded = simulated band)')
ax[0,1].set_title('loss tangent')

for n, (name, fit) in enumerate(FITS.items()):
    r = results[name]
    e = eps_poles(freq, fit)
    ax[1,0].plot(freq/1e9, r['eps_eff'], f'C{n}-', label=name)
    ax[1,0].plot(freq/1e9, hammerstad(e.real), f'C{n}--', lw=1)
    ax[1,1].plot(freq/1e9, r['alpha_dB'], f'C{n}-', label=name)
    ax[1,1].plot(freq/1e9, alpha_d_dB_m(e.real, r['eps_eff'], -e.imag/e.real, freq),
                 f'C{n}--', lw=1)
ax[1,0].set_ylabel('eps_eff of the line'); ax[1,1].set_ylabel('alpha (dB/m)')
ax[1,0].set_title('solid: openEMS   dashed: closed form')
ax[1,1].set_title('dielectric attenuation')
for a in (ax[1,0], ax[1,1]):
    a.set_xlabel('f (GHz)'); a.grid(); a.legend(fontsize=8)

plt.show()
