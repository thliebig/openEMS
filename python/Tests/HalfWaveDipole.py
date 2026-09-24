# -*- coding: utf-8 -*-
"""
 Half-wave dipole — CurvePort feed test

 Verifies CurvePort by feeding a λ/2 wire dipole (arm = 75 mm, total ≈ 150 mm)
 at its center gap.  Checks that the FDTD result reproduces known thin-wire
 analytic values:

 Pass criteria:
   resonant frequency within 10 % of 0.94 × c/(2L)  (≈ 0.94 GHz)
   Re(Z_in) at resonance within 20 % of 73 Ω
   min(S11) < −15 dB at resonance  (73 Ω reference)

 Tested with
  - python 3.14
  - openEMS v0.0.36+

 (c) 2026 Thorsten Liebig <thorsten.liebig@gmx.de>

"""

import os, tempfile
import numpy as np

from CSXCAD  import ContinuousStructure
from openEMS import openEMS
from openEMS.physical_constants import C0
from openEMS.ports import CurvePort

### Simulation path
Sim_Path = os.path.join(tempfile.gettempdir(), 'HalfWaveDipole')

### Geometry
unit     = 1e-3   # drawing unit: mm
arm_len  = 75     # each dipole arm (mm); total ≈ 150 mm
mesh_res = 5.0    # mm

### Frequency
f_start   = 0.7e9
f_stop    = 1.3e9
f_0       = 1.0e9
lambda0   = C0 / f_0 / unit   # 300 mm free-space wavelength at f_0

# A dipole first resonates at L ≈ 0.47λ (not λ/2); the 0.94 factor corrects
# for end effects and feed-gap capacitance loading.
f_res_est = C0 / (2 * 2*arm_len * unit) * 0.94   # ≈ 0.94 GHz

### FDTD
FDTD = openEMS(EndCriteria=1e-4)
FDTD.SetGaussExcite(0.5*(f_start + f_stop), 0.5*(f_stop - f_start))
FDTD.SetBoundaryCond(['PML_8'] * 6)

### CSXCAD geometry & mesh
CSX = ContinuousStructure()
FDTD.SetCSX(CSX)
mesh = CSX.GetGrid()
mesh.SetDeltaUnit(unit)

pad = lambda0 / 4   # 75 mm — clearance between dipole tips and PML
pad += 10*mesh_res  # add space for the PML

mesh.AddLine('x', [-pad, 0, pad])
mesh.AddLine('y', [-pad, 0, pad])
mesh.AddLine('z', [-(arm_len + pad), 0, arm_len + pad])
mesh.SmoothMeshLines('all', mesh_res, ratio=1.4)

### CurvePort spanning the full dipole length.
### Internally it places the single-cell feed at the midpoint (z = 0) and
### adds PEC arm wires from each endpoint to the feed cell automatically.
port = CurvePort(CSX, 1, R=73, start=[0, 0, -arm_len], stop=[0, 0, arm_len], excite=1)

if 0:  # debugging only
    CSX_file = os.path.join(Sim_Path, 'dipole.xml')
    if not os.path.exists(Sim_Path):
        os.mkdir(Sim_Path)
    CSX.Write2XML(CSX_file)
    from CSXCAD import AppCSXCAD_BIN
    os.system(AppCSXCAD_BIN + ' "{}"'.format(CSX_file))

### Run
FDTD.Run(Sim_Path, cleanup=True, exact_endcriteria=True)

### Post-processing
freq = np.linspace(f_start, f_stop, 401)
port.CalcPort(Sim_Path, freq)

s11    = port.uf_ref / port.uf_inc
s11_dB = 20 * np.log10(np.abs(s11))
Z_in   = port.Z_ref * (1 + s11) / (1 - s11)

# Resonance = S11 minimum in band
i_res = np.argmin(s11_dB)
f_res = freq[i_res]
Z_res = np.real(Z_in[i_res])

### Pass / fail checks
print(f'Resonant frequency:    {f_res/1e9:.3f} GHz  (corrected estimate: {f_res_est/1e9:.3f} GHz)')
print(f'Re(Z_in) at resonance: {Z_res:.1f} Ω  (analytic: ~73 Ω)')
print(f'min(S11) = {s11_dB[i_res]:.1f} dB')

assert f_start < f_res < f_stop, \
    f'FAIL: resonance at {f_res/1e9:.3f} GHz is outside the sweep ' \
    f'[{f_start/1e9:.1f}, {f_stop/1e9:.1f}] GHz'
assert abs(f_res - f_res_est) / f_res_est < 0.10, \
    f'FAIL: f_res = {f_res/1e9:.3f} GHz, resonance estimate = {f_res_est/1e9:.3f} GHz (>10 % off)'
assert abs(Z_res - 73) / 73 < 0.20, \
    f'FAIL: Re(Z_in) = {Z_res:.1f} Ω, expected 73 Ω (±20 %)'
assert s11_dB[i_res] < -15, \
    f'FAIL: min(S11) = {s11_dB[i_res]:.1f} dB, expected < -15 dB at resonance'

print('PASS')

if 0:  # set to 1 for debugging plots
    import matplotlib.pyplot as plt

    fig, axis = plt.subplots(num='S11', tight_layout=True)
    axis.plot(freq/1e9, s11_dB, 'k-', linewidth=2, label='$S_{11}$')
    axis.axvline(f_res/1e9,     color='r', linestyle='--', label=f'f_res = {f_res/1e9:.3f} GHz')
    axis.axvline(f_res_est/1e9, color='g', linestyle=':',  label=f'f_est = {f_res_est/1e9:.3f} GHz')
    axis.grid()
    axis.set_xmargin(0)
    axis.set_xlabel('Frequency (GHz)')
    axis.set_ylabel('$S_{11}$ (dB)')
    axis.legend()

    fig, axes = plt.subplots(2, 1, num='Input Impedance', tight_layout=True, sharex=True)
    axes[0].plot(freq/1e9, np.real(Z_in), 'k-',  linewidth=2, label=r'$\Re\{Z_{in}\}$')
    axes[0].axhline(73, color='g', linestyle=':', linewidth=1, label='73 Ω (analytic)')
    axes[0].set_ylabel(r'$\Re\{Z_{in}\}$ (Ω)')
    axes[0].legend()
    axes[0].grid()
    axes[1].plot(freq/1e9, np.imag(Z_in), 'r--', linewidth=2, label=r'$\Im\{Z_{in}\}$')
    axes[1].axhline(0, color='gray', linestyle=':', linewidth=1)
    axes[1].set_ylabel(r'$\Im\{Z_{in}\}$ (Ω)')
    axes[1].set_xlabel('Frequency (GHz)')
    axes[1].legend()
    axes[1].grid()
    axes[1].set_xmargin(0)

    plt.show()
