# -*- coding: utf-8 -*-
"""
 Stripline Test (PML-terminated)

 Verifies StripLinePort setup and S-parameter extraction for a stripline
 embedded in a substrate between two PEC ground planes (enforced by the
 top/bottom PEC boundary conditions).  PML absorbers terminate the x-direction.

 Pass criteria:
   max(dB(S11)) < -40 dB   (low reflection)
   min(dB(S21)) > -0.1 dB  (near-lossless transmission)

 Tested with
  - python 3.14
  - openEMS v0.0.36+

 (c) 2026 Thorsten Liebig <thorsten.liebig@gmx.de>

"""

import os, tempfile
import numpy as np

from CSXCAD  import ContinuousStructure
from openEMS import openEMS
from openEMS.physical_constants import *
from openEMS.ports import StripLinePort

### Setup the simulation
Sim_Path = os.path.join(tempfile.gettempdir(), 'Stripline')

unit          = 1e-6  # drawing unit in um
SL_length     = 50000
SL_width      = 520
SL_height     = 500   # distance from strip to each ground plane
substrate_epr = 3.66
f_max         = 7e9

### Setup FDTD parameters & excitation
FDTD = openEMS(EndCriteria=1e-4)
FDTD.SetGaussExcite(f_max/2, f_max/2)
FDTD.SetBoundaryCond(['PML_8', 'PML_8', 'PMC', 'PMC', 'PEC', 'PEC'])

### Setup CSXCAD geometry & mesh
CSX = ContinuousStructure()
FDTD.SetCSX(CSX)
mesh = CSX.GetGrid()
mesh.SetDeltaUnit(unit)

resolution = C0 / (f_max * np.sqrt(substrate_epr)) / unit / 50
third_mesh = np.array([-1/3, 2/3]) * resolution / 4

mesh.AddLine('x', [-SL_length/2, 0, SL_length/2])
mesh.SmoothMeshLines('x', resolution, ratio=1.5)

mesh.AddLine('y', [0])
mesh.AddLine('y', SL_width/2 + third_mesh)
mesh.SmoothMeshLines('y', resolution/4, ratio=1.5)
y_pos = mesh.GetLines('y')
mesh.AddLine('y', np.concatenate([-y_pos, [-10*SL_width, 10*SL_width]]))
mesh.SmoothMeshLines('y', resolution, ratio=1.3)

mesh.AddLine('z', np.linspace(0, SL_height, 5))
mesh.AddLine('z', -np.linspace(0, SL_height, 5))

### Substrate
substrate = CSX.AddMaterial('RO4350B', epsilon=substrate_epr)
start = [mesh.GetLines('x')[0],  mesh.GetLines('y')[0],  mesh.GetLines('z')[0]]
stop  = [mesh.GetLines('x')[-1], mesh.GetLines('y')[-1], mesh.GetLines('z')[-1]]
substrate.AddBox(start, stop)

### Stripline ports (include the metal strip)
pec = CSX.AddMetal('PEC')

portstart = [mesh.GetLines('x')[0], -SL_width/2, 0]
portstop  = [0,                      SL_width/2, 0]
port1 = StripLinePort(CSX, 1, pec, portstart, portstop, 'x', 'z', SL_height,
                      excite=1, priority=999,
                      FeedShift=10*resolution, MeasPlaneShift=SL_length/3)

portstart = [mesh.GetLines('x')[-1], -SL_width/2, 0]
portstop  = [0,                       SL_width/2, 0]
port2 = StripLinePort(CSX, 2, pec, portstart, portstop, 'x', 'z', SL_height,
                      priority=999, MeasPlaneShift=SL_length/3)

ports = [port1, port2]

### Run the simulation
FDTD.Run(Sim_Path, cleanup=True, exact_endcriteria=True)

### Post-processing
f = np.linspace(1e6, f_max, 1601)
for port in ports:
    port.CalcPort(Sim_Path, f, ref_impedance=50)

s11 = ports[0].uf_ref / ports[0].uf_inc
s21 = ports[1].uf_ref / ports[0].uf_inc

s11_dB = 20 * np.log10(np.abs(s11))
s21_dB = 20 * np.log10(np.abs(s21))

### Pass / fail checks
mask = f > 100e6
print(f'max(dB(S11)) = {np.max(s11_dB[mask]):.1f} dB')
print(f'min(dB(S21)) = {np.min(s21_dB[mask]):.1f} dB,  max(dB(S21)) = {np.max(s21_dB[mask]):.2f} dB')

assert np.max(s11_dB[mask]) < -40, \
    f'FAIL: max(dB(S11)) = {np.max(s11_dB[mask]):.1f} dB, expected < -40 dB'
assert np.min(s21_dB[mask]) > -0.1, \
    f'FAIL: min(dB(S21)) = {np.min(s21_dB[mask]):.1f} dB, expected > -0.1 dB'
assert np.max(s21_dB[mask]) < 0.05, \
    f'FAIL: max(dB(S21)) = {np.max(s21_dB[mask]):.2f} dB, expected < +0.05 dB (sign error?)'

print('PASS')

if 0:  # set to 1 for debugging plots
    import matplotlib.pyplot as plt

    fig, axis = plt.subplots(num='S-Parameters', tight_layout=True)
    axis.plot(f/1e9, s11_dB, 'k-',  linewidth=2, label='$S_{11}$')
    axis.plot(f/1e9, s21_dB, 'r--', linewidth=2, label='$S_{21}$')
    axis.grid()
    axis.set_xmargin(0)
    axis.set_ylim([-50, 2])
    axis.set_xlabel('Frequency (GHz)')
    axis.set_ylabel('S-Parameter (dB)')
    axis.legend()

    plt.show()
