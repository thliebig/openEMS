# -*- coding: utf-8 -*-
"""
 Conducting sheet test — plane-wave transmission through a thin metal sheet

 A TEM plane wave (Ex, travelling in +z) is launched in a 1D channel built
 from PEC walls normal to x, PMC walls normal to y and PML at both z-ends.
 A conducting sheet (zero-thickness model of a metal layer with finite
 conductivity sigma and thickness t) spans the full channel cross-section.

 The conducting sheet model approximates the sheet admittance

     Y(f) = 2 sigma / gamma * tanh(gamma t / 2),   gamma = sqrt(j w mue0 sigma)

 so the transmission, normalised to a reference run without the sheet, is

     T(f) = 2 / (2 + Z0 Y(f))

 Test cases
 ----------
 1. Resistive film   sigma = 1e5 S/m,   t = 10 um  (t < skin depth, flat response)
 2. Thin copper      sigma = 5.8e7 S/m, t = 1 um   (t ~ skin depth at 10 GHz)
 3. Thick copper     sigma = 5.8e7 S/m, t = 18 um  (t >> skin depth, T ~ sqrt(f))

 Pass criteria (1 .. 10 GHz)
   case 1, 2 : |dB(T_fdtd) - dB(T_analytic)| < 0.25 dB,  |arg(T_fdtd / T_analytic)| < 10 deg
   case 3    : |dB(T_fdtd) - dB(T_analytic)| < 2.5 dB  (the model accuracy degrades
               for t >> skin depth, this guards against gross errors only)
   case 3    : T rises with frequency by 10 dB/decade (skin effect) within 2 dB

 Tested with
  - python 3.13
  - openEMS v0.0.37+

 (c) 2026 Sean Mollet <sean@malmoset.com>

"""

import os, tempfile
import numpy as np

from CSXCAD  import ContinuousStructure
from openEMS import openEMS
from openEMS.physical_constants import Z0, MUE0
from openEMS.utilities import DFT_time2freq

### Geometry
unit     = 1e-3    # drawing unit: mm
dz       = 0.25    # axial mesh resolution (mm)
width    = 1.0     # transverse channel width (mm)
length   = 150     # channel length (mm)
z_src    = 20      # excitation plane
z_sheet  = 60      # conducting sheet position
z_probe  = 110     # transmitted-field probe

### Frequency
f_start = 1e9
f_stop  = 10e9
freq    = np.linspace(f_start, f_stop, 91)


def run(Sim_Path, sheet=None):
    """ Run the channel with an optional (conductivity, thickness) sheet and return Ex(freq) at the probe. """
    FDTD = openEMS(NrTS=20000, EndCriteria=1e-5)
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

    if sheet is not None:
        cs = CSX.AddConductingSheet('sheet', conductivity=sheet[0], thickness=sheet[1])
        cs.AddBox([0, 0, z_sheet], [width, width, z_sheet], priority=10)

    FDTD.Run(Sim_Path, cleanup=True)

    # field probe file columns: t/s, Ex, Ey, Ez
    data = np.loadtxt(os.path.join(Sim_Path, 'et_trans'), comments='%')
    return DFT_time2freq(data[:, 0], data[:, 1], freq)


def sheet_transmission(sigma, thickness):
    gamma = np.sqrt(1j * 2*np.pi*freq * MUE0 * sigma)
    Y = 2*sigma/gamma * np.tanh(gamma*thickness/2)
    return 2 / (2 + Z0*Y)


### Reference run (empty channel)
print('Running reference (no sheet)...')
E_ref = run(os.path.join(tempfile.gettempdir(), 'CondSheet_Ref'))

cases = [('Resistive film', 1e5,   10e-6, 0.25, 10),
         ('Thin copper',    5.8e7,  1e-6, 0.25, 10),
         ('Thick copper',   5.8e7, 18e-6, 2.5,  None)]

results = {}
for name, sigma, thickness, dB_tol, phase_tol in cases:
    print('Running {} (sigma = {:g} S/m, t = {:g} um)...'.format(name, sigma, thickness*1e6))
    Sim_Path = os.path.join(tempfile.gettempdir(), 'CondSheet_' + name.replace(' ', '_'))

    T_fdtd = run(Sim_Path, (sigma, thickness)) / E_ref
    T_ana  = sheet_transmission(sigma, thickness)
    results[name] = (T_fdtd, T_ana)

    dB_err    = 20*np.log10(np.abs(T_fdtd / T_ana))
    phase_err = np.angle(T_fdtd / T_ana, deg=True)
    print(f'  |T| = {20*np.log10(np.abs(T_fdtd[0])):.1f} .. {20*np.log10(np.abs(T_fdtd[-1])):.1f} dB'
          f'  (analytic {20*np.log10(np.abs(T_ana[0])):.1f} .. {20*np.log10(np.abs(T_ana[-1])):.1f} dB)')
    print(f'  max |dB error| = {np.max(np.abs(dB_err)):.3f} dB,  max |phase error| = {np.max(np.abs(phase_err)):.2f} deg')

    assert np.max(np.abs(dB_err)) < dB_tol, \
        f'FAIL [{name}]: |T| deviates {np.max(np.abs(dB_err)):.3f} dB from analytic, expected < {dB_tol} dB'
    if phase_tol is not None:
        assert np.max(np.abs(phase_err)) < phase_tol, \
            f'FAIL [{name}]: arg(T) deviates {np.max(np.abs(phase_err)):.2f} deg from analytic, expected < {phase_tol} deg'

    print('PASS [{}]'.format(name))

### Skin effect: for t >> skin depth |Y| ~ 1/sqrt(f), i.e. |T| rises by 10 dB/decade
T_fdtd, _ = results['Thick copper']
slope = 20*np.log10(np.abs(T_fdtd[-1] / T_fdtd[0])) / np.log10(freq[-1] / freq[0])
print(f'Thick copper |T| slope: {slope:.1f} dB/decade (skin effect: 10 dB/decade)')
assert abs(slope - 10) < 2, \
    f'FAIL [Thick copper]: |T| slope {slope:.1f} dB/decade, expected 10 +/- 2 dB/decade'

print('PASS')

if 0:  # set to 1 for debugging plots
    import matplotlib.pyplot as plt

    fig, axis = plt.subplots(num='Sheet transmission', tight_layout=True)
    for n, (name, (T_fdtd, T_ana)) in enumerate(results.items()):
        axis.plot(freq/1e9, 20*np.log10(np.abs(T_fdtd)), f'C{n}-',  linewidth=2, label=name)
        axis.plot(freq/1e9, 20*np.log10(np.abs(T_ana)),  f'C{n}--', linewidth=1)
    axis.grid()
    axis.set_xmargin(0)
    axis.set_xlabel('Frequency (GHz)')
    axis.set_ylabel('|T| (dB)')
    axis.legend()
    plt.show()
