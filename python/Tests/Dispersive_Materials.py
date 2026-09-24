# -*- coding: utf-8 -*-
"""
 Dispersive material test — plane-wave transmission through a slab

 A TEM plane wave (Ex, travelling in +z) is launched in a 1D channel built
 from PEC walls normal to x, PMC walls normal to y and PML at both z-ends.
 A slab of dispersive material fills the full channel cross-section.  The
 transmitted field, normalised to a reference run without the slab, is
 compared to the analytic slab transmission

     T(f) = exp(j k0 d) / (cos(k d) + j/2 (eta + 1/eta) sin(k d))

 with k = k0 sqrt(eps_r mue_r) and eta = sqrt(mue_r / eps_r)  (e^{jwt} convention).

 Test cases
 ----------
 1. Drude       (electric)                  eps_r(f) = 1 - fp^2/(f^2 - j f/(2 pi tau))
 2. Lorentz     (electric, pole at f_pole)  eps_r(f) = 1 - fp^2/(f^2 - f_pole^2 - j f/(2 pi tau))
 3. Debye       (electric)                  eps_r(f) = eps_inf + d_eps/(1 + j w tau)
 4. Double Drude (electric + magnetic)      eps_r = mue_r --> eta = 1, reflection free,
                                            negative refractive index below fp/sqrt(2),
                                            nearly lossless

 Pass criteria (per case, 2 .. 9 GHz)
   max |T_fdtd - T_analytic| < 0.05
   Double Drude: max(dB(|T|)) deviation < 0.5 dB (no reflection at any frequency)

 Tested with
  - python 3.13
  - openEMS v0.0.37+

 (c) 2026 Sean Mollet <sean@malmoset.com>

"""

import os, tempfile
import numpy as np

from CSXCAD  import ContinuousStructure
from CSXCAD.CSProperties import CSPropLorentzMaterial, CSPropDebyeMaterial
from openEMS import openEMS
from openEMS.physical_constants import C0
from openEMS.utilities import DFT_time2freq

### Geometry
unit     = 1e-3    # drawing unit: mm
dz       = 0.25    # axial mesh resolution (mm)
width    = 1.0     # transverse channel width (mm)
length   = 150     # channel length (mm)
z_src    = 20      # excitation plane
slab_z0  = 60      # slab start
slab_d   = 10      # slab thickness (mm)
z_probe  = 110     # transmitted-field probe

### Frequency
f_start = 1e9
f_stop  = 10e9
freq    = np.linspace(2e9, 9e9, 141)   # evaluation band
w       = 2*np.pi*freq


def run(Sim_Path, material_fn=None):
    """ Run the channel with an optional slab and return the probe spectrum Ex(freq). """
    FDTD = openEMS(NrTS=20000, EndCriteria=1e-6)
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

    if material_fn is not None:
        mat = material_fn(CSX)
        mat.AddBox([0, 0, slab_z0], [width, width, slab_z0 + slab_d], priority=10)

    FDTD.Run(Sim_Path, cleanup=True)

    # field probe file columns: t/s, Ex, Ey, Ez
    data = np.loadtxt(os.path.join(Sim_Path, 'et_trans'), comments='%')
    return DFT_time2freq(data[:, 0], data[:, 1], freq)


def slab_transmission(eps_r, mue_r):
    k0  = w / C0
    k   = k0 * np.sqrt(eps_r * mue_r + 0j)
    k   = np.where(np.imag(k) > 0, -k, k)   # passive medium: Im(k) <= 0 for e^{jwt}
    eta = np.sqrt(mue_r / eps_r + 0j)
    eta = np.where(np.real(eta) < 0, -eta, eta)
    d   = slab_d * unit
    return np.exp(1j*k0*d) / (np.cos(k*d) + 0.5j*(eta + 1/eta)*np.sin(k*d))


### Reference run (empty channel)
print('Running reference (no slab)...')
E_ref = run(os.path.join(tempfile.gettempdir(), 'Dispersive_Ref'))


### Material definitions
fp   = 5e9        # plasma frequency (Hz)
tau  = 1e-9       # relaxation time (s)

def drude(CSX):
    m = CSPropLorentzMaterial(CSX.GetParameterSet(), order=1)
    m.SetName('drude')
    m.SetDispersiveMaterialProperty(0, eps_plasma=fp, eps_relax=tau)
    CSX.AddProperty(m)
    return m

def eps_drude(f, tau=tau):
    return 1 - fp**2 / (f**2 - 1j*f/(2*np.pi*tau))

f_pole = 3e9
def lorentz(CSX):
    m = CSPropLorentzMaterial(CSX.GetParameterSet(), order=1)
    m.SetName('lorentz')
    m.SetDispersiveMaterialProperty(0, eps_plasma=fp, eps_pole_freq=f_pole, eps_relax=tau)
    CSX.AddProperty(m)
    return m

def eps_lorentz(f):
    return 1 - fp**2 / (f**2 - f_pole**2 - 1j*f/(2*np.pi*tau))

# Note: the Debye ADE becomes unstable for d_eps >~ eps_inf, keep d_eps small
eps_inf   = 4.0
d_eps     = 1.0
tau_debye = 1 / (2*np.pi*4e9)   # relaxation frequency 4 GHz
def debye(CSX):
    m = CSPropDebyeMaterial(CSX.GetParameterSet(), order=1, epsilon=eps_inf)
    m.SetName('debye')
    m.SetDispersiveMaterialProperty(0, eps_delta=d_eps, eps_relax=tau_debye)
    CSX.AddProperty(m)
    return m

def eps_debye(f):
    return eps_inf + d_eps / (1 + 1j*2*np.pi*f*tau_debye)

tau_mtm = 1e-7    # nearly lossless, |T| stays within 0.02 dB of 0 dB
def double_drude(CSX):
    m = CSPropLorentzMaterial(CSX.GetParameterSet(), order=1)
    m.SetName('double_drude')
    m.SetDispersiveMaterialProperty(0, eps_plasma=fp, eps_relax=tau_mtm,
                                       mue_plasma=fp, mue_relax=tau_mtm)
    CSX.AddProperty(m)
    return m


cases = [('Drude',        drude,        eps_drude(freq),   np.ones_like(freq)),
         ('Lorentz',      lorentz,      eps_lorentz(freq), np.ones_like(freq)),
         ('Debye',        debye,        eps_debye(freq),   np.ones_like(freq)),
         ('Double Drude', double_drude, eps_drude(freq, tau_mtm), eps_drude(freq, tau_mtm))]

results = {}
for name, material_fn, eps_r, mue_r in cases:
    print('Running {}...'.format(name))
    Sim_Path = os.path.join(tempfile.gettempdir(), 'Dispersive_' + name.replace(' ', '_'))
    E_slab = run(Sim_Path, material_fn)

    T_fdtd = E_slab / E_ref
    T_ana  = slab_transmission(eps_r, mue_r)
    results[name] = (T_fdtd, T_ana)

    err = np.max(np.abs(T_fdtd - T_ana))
    print(f'  max |T_fdtd - T_analytic| = {err:.4f}')
    assert err < 0.05, \
        f'FAIL [{name}]: max |T_fdtd - T_analytic| = {err:.4f}, expected < 0.05'

    if name == 'Double Drude':
        T_dB = 20*np.log10(np.abs(T_fdtd))
        print(f'  |T| range: {np.min(T_dB):.2f} .. {np.max(T_dB):.2f} dB')
        assert np.max(np.abs(T_dB)) < 0.5, \
            f'FAIL [{name}]: |T| deviates {np.max(np.abs(T_dB)):.2f} dB from 0 dB, expected < 0.5 dB (eta=1 slab must not reflect)'

    print('PASS [{}]'.format(name))

print('PASS')

if 0:  # set to 1 for debugging plots
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 1, num='Slab transmission', tight_layout=True, sharex=True)
    for n, (name, (T_fdtd, T_ana)) in enumerate(results.items()):
        axes[0].plot(freq/1e9, 20*np.log10(np.abs(T_fdtd)), f'C{n}-',  linewidth=2, label=name)
        axes[0].plot(freq/1e9, 20*np.log10(np.abs(T_ana)),  f'C{n}--', linewidth=1)
        axes[1].plot(freq/1e9, np.angle(T_fdtd, deg=True),  f'C{n}-',  linewidth=2)
        axes[1].plot(freq/1e9, np.angle(T_ana,  deg=True),  f'C{n}--', linewidth=1)
    axes[0].set_ylabel('|T| (dB)')
    axes[1].set_ylabel('arg(T) (deg)')
    axes[1].set_xlabel('Frequency (GHz)')
    for ax in axes:
        ax.grid()
        ax.set_xmargin(0)
    axes[0].legend()
    plt.show()
