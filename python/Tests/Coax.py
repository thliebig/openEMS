# -*- coding: utf-8 -*-
"""
 Coaxial Line Test

 Verifies CoaxialPort setup, S-parameter extraction, and line impedance
 against the analytic formula for a 50 Ω air-filled coaxial line.

 Pass criteria:
   max(dB(S11)) < -0.1 dB   (adequate return loss over band)
   min(dB(S21)) > -50 dB    (cable transmits signal)
   Z_ref within 5 % of analytic value (mid-band)

 Tested with
  - python 3.10
  - openEMS v0.0.35+

 (c) 2026 Thorsten Liebig <thorsten.liebig@gmx.de>

"""

import os, tempfile
import numpy as np

from CSXCAD  import ContinuousStructure
from openEMS import openEMS
from openEMS.physical_constants import *
from openEMS.ports import CoaxialPort

### Setup the simulation
Sim_Path = os.path.join(tempfile.gettempdir(), 'Coax')

unit = 1e-3          # drawing unit in mm
length       = 150   # cable length (mm)
coax_rad_i   = 5     # inner conductor radius  (5 mm)
coax_rad_ai  = 11.5  # inner radius of outer conductor  (ratio 2.3 → 50 Ω)
coax_rad_aa  = 12    # outer radius of outer conductor  (12 mm)
mesh_res     = [0.25, 0.25, 2]  # 0.25 mm transverse, 2 mm axial

f_stop = 3e9
epsR   = 1

### Setup FDTD parameters & excitation
FDTD = openEMS(EndCriteria=1e-4)
FDTD.SetGaussExcite(0, f_stop)
FDTD.SetBoundaryCond(['PEC', 'PEC', 'PEC', 'PEC', 'PML_8', 'PML_8'])

### Setup CSXCAD geometry & mesh
CSX = ContinuousStructure()
FDTD.SetCSX(CSX)
mesh = CSX.GetGrid()
mesh.SetDeltaUnit(unit)

mesh.AddLine('x', np.arange(-coax_rad_aa, coax_rad_aa + mesh_res[0], mesh_res[0]))
mesh.AddLine('y', np.arange(-coax_rad_aa, coax_rad_aa + mesh_res[1], mesh_res[1]))
mesh.AddLine('z', np.arange(0, length + mesh_res[2],                 mesh_res[2]))

### Coaxial ports
copper = CSX.AddMetal('copper')
ports  = []

start = [0, 0, 0]
stop  = [0, 0, length/2]
ports.append(CoaxialPort(CSX, 1, copper, None, start, stop, 'z',
                         coax_rad_i, coax_rad_ai, coax_rad_aa,
                         excite=1, FeedShift=10*mesh_res[2], priority=10))

start = [0, 0, length]
stop  = [0, 0, length/2]
ports.append(CoaxialPort(CSX, 2, copper, None, start, stop, 'z',
                         coax_rad_i, coax_rad_ai, coax_rad_aa, priority=10))

### Run the simulation
FDTD.Run(Sim_Path, cleanup=True, exact_endcriteria=True)

### Post-processing
freq = np.linspace(1e6, f_stop, 201)
for port in ports:
    port.CalcPort(Sim_Path, freq)

s11 = ports[0].uf_ref / ports[0].uf_inc
s21 = ports[1].uf_ref / ports[0].uf_inc

s11_dB = 20 * np.log10(np.abs(s11))
s21_dB = 20 * np.log10(np.abs(s21))

### Pass / fail checks (skip the DC-vicinity bin)
mask = freq > 100e6
print(f'max(dB(S11)) = {np.max(s11_dB[mask]):.1f} dB')
print(f'min(dB(S21)) = {np.min(s21_dB[mask]):.1f} dB,  max(dB(S21)) = {np.max(s21_dB[mask]):.2f} dB')

assert np.max(s11_dB[mask]) < -50, \
    f'FAIL: max(dB(S11)) = {np.max(s11_dB[mask]):.1f} dB, expected < -50 dB'
assert np.min(s21_dB[mask]) > -0.1, \
    f'FAIL: min(dB(S21)) = {np.min(s21_dB[mask]):.1f} dB, expected > -0.1 dB'
assert np.max(s21_dB[mask]) < 0.05, \
    f'FAIL: max(dB(S21)) = {np.max(s21_dB[mask]):.2f} dB, expected < +0.05 dB (sign error?)'

ZL_a   = Z0 / (2*np.pi) / np.sqrt(epsR) * np.log(coax_rad_ai / coax_rad_i)
ZL_num = np.real(ports[0].Z_ref[mask])
rel_err = np.max(np.abs(ZL_num - ZL_a) / ZL_a)
assert rel_err < 0.05, \
    f'FAIL: line impedance error {rel_err*100:.1f} %, expected < 5 %'

print('PASS')

if 0:  # set to 1 for debugging plots
    import matplotlib.pyplot as plt

    fig, axis = plt.subplots(num='S-Parameters', tight_layout=True)
    axis.plot(freq/1e9, s11_dB, linewidth=2, label='$S_{11}$')
    axis.plot(freq/1e9, s21_dB, 'r--', linewidth=2, label='$S_{21}$')
    axis.grid()
    axis.set_xmargin(0)
    axis.set_xlabel('Frequency (GHz)')
    axis.set_ylabel('S-Parameter (dB)')
    axis.legend()

    fig, axis = plt.subplots(num='Line Impedance', tight_layout=True)
    axis.plot(freq/1e9, np.real(ports[0].Z_ref), linewidth=2, label=r'$\Re\{Z_L\}$')
    axis.plot(freq/1e9, np.imag(ports[0].Z_ref), 'r--', linewidth=2, label=r'$\Im\{Z_L\}$')
    axis.plot(freq/1e9, ZL_a * np.ones_like(freq), 'g-.', linewidth=2, label='$Z_{L,analytic}$')
    axis.grid()
    axis.set_xmargin(0)
    axis.set_xlabel('Frequency (GHz)')
    axis.set_ylabel(r'Line impedance $(\Omega)$')
    axis.legend()

    plt.show()
