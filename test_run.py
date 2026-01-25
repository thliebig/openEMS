#!/usr/bin/env python3
"""
Simple openEMS test - Rectangular Waveguide
"""

import os, tempfile
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-GUI backend
from matplotlib import pyplot as plt

from CSXCAD import ContinuousStructure
from openEMS import openEMS
from openEMS.physical_constants import C0

# Setup the simulation
Sim_Path = os.path.join(tempfile.gettempdir(), 'Rect_WG_Test')

unit = 1e-6  # drawing unit in um

# waveguide dimensions (WR42)
a = 10700    # waveguide width
b = 4300     # waveguide height
length = 50000

# frequency range of interest
f_start = 20e9
f_0     = 24e9
f_stop  = 26e9
lambda0 = C0/f_0/unit

# waveguide TE-mode definition
TE_mode = 'TE10'

# targeted mesh resolution
mesh_res = lambda0/30

### Setup FDTD parameter & excitation function
FDTD = openEMS(NrTS=1e4)
FDTD.SetGaussExcite(0.5*(f_start+f_stop), 0.5*(f_stop-f_start))

# boundary conditions
FDTD.SetBoundaryCond([0, 0, 0, 0, 3, 3])

### Setup geometry & mesh
CSX = ContinuousStructure()
FDTD.SetCSX(CSX)
mesh = CSX.GetGrid()
mesh.SetDeltaUnit(unit)

mesh.AddLine('x', [0, a])
mesh.AddLine('y', [0, b])
mesh.AddLine('z', [0, length])

## Apply the waveguide port
ports = []
start = [0, 0, 10*mesh_res]
stop  = [a, b, 15*mesh_res]
mesh.AddLine('z', [start[2], stop[2]])
ports.append(FDTD.AddRectWaveGuidePort(0, start, stop, 'z', a*unit, b*unit, TE_mode, 1))

start = [0, 0, length-10*mesh_res]
stop  = [a, b, length-15*mesh_res]
mesh.AddLine('z', [start[2], stop[2]])
ports.append(FDTD.AddRectWaveGuidePort(1, start, stop, 'z', a*unit, b*unit, TE_mode))

mesh.SmoothMeshLines('all', mesh_res, ratio=1.4)

### Define dump box
Et = CSX.AddDump('Et', file_type=0, sub_sampling=[2,2,2])
start = [0, 0, 0]
stop  = [a, b, length]
Et.AddBox(start, stop)

### Run the simulation
print("Starting openEMS simulation...")
FDTD.Run(Sim_Path, cleanup=True)
print("Simulation completed!")

### Postprocessing
freq = np.linspace(f_start, f_stop, 201)
for port in ports:
    port.CalcPort(Sim_Path, freq)

s11 = ports[0].uf_ref / ports[0].uf_inc
s21 = ports[1].uf_ref / ports[0].uf_inc

## Plot s-parameter
plt.figure()
plt.plot(freq*1e-6, 20*np.log10(np.abs(s11)), 'k-', linewidth=2, label='$S_{11}$')
plt.grid()
plt.plot(freq*1e-6, 20*np.log10(np.abs(s21)), 'r--', linewidth=2, label='$S_{21}$')
plt.legend()
plt.ylabel('S-Parameter (dB)')
plt.xlabel('frequency (MHz)')
plt.savefig('/home/user/openems/test_results.png')
print("Results saved to test_results.png")

# Print some results
print(f"\nS11 at {f_0/1e9:.1f} GHz: {20*np.log10(np.abs(s11[100])):.2f} dB")
print(f"S21 at {f_0/1e9:.1f} GHz: {20*np.log10(np.abs(s21[100])):.2f} dB")
print("\nTest completed successfully!")
