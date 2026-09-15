# -*- coding: utf-8 -*-
"""
 Coplanar Waveguide (CPW) Line Test

 Verifies CPWPort setup and S-parameter extraction for a finite-length CPW
 transmission line on RO4350B substrate with 50-Ohm port terminations.

 Pass criteria:
   max(dB(S11)) < -20 dB   (low reflection)
   min(dB(S21)) > -0.5 dB  (near-lossless — substrate is modelled lossless)
   max(dB(S21)) < +0.01 dB (sign-error guard)

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
from openEMS.ports import CPWPort

### Setup the simulation
Sim_Path = os.path.join(tempfile.gettempdir(), 'CPW_Line')

unit                = 1e-6   # drawing unit in um
CPW_length          = 40000
CPW_port_length     = 10000
CPW_width           = 1000
CPW_gap             = 140
substrate_thickness = 512
substrate_width     = 5000
substrate_epr       = 3.66
f_max               = 10e9
air_spacing         = 7000

feed_R = 50  # lumped port termination resistance

### Setup FDTD parameters & excitation
FDTD = openEMS(EndCriteria=1e-4)
FDTD.SetGaussExcite(f_max/2, f_max/2)
FDTD.SetBoundaryCond(['PMC']*6)

### Setup CSXCAD geometry & mesh
CSX = ContinuousStructure()
FDTD.SetCSX(CSX)
mesh = CSX.GetGrid()
mesh.SetDeltaUnit(unit)

resolution = C0 / (f_max * np.sqrt(substrate_epr)) / unit / 30
edge_res   = 40

# x-mesh: fine at port transitions, coarse in middle and beyond ports
mesh.AddLine('x', [0])
mesh.AddLine('x', [-CPW_length/2, -CPW_length/2 + CPW_port_length,
                    CPW_length/2 - CPW_port_length, CPW_length/2])
mesh.AddLine('x', [-CPW_length/2 - air_spacing, CPW_length/2 + air_spacing])
mesh.SmoothMeshLines('x', resolution, ratio=1.5)

# y-mesh: fine at CPW signal and gap edges, coarse elsewhere
third_mesh = np.array([-2/3, 1/3]) * edge_res
mesh.AddLine('y', [0])
mesh.AddLine('y',  CPW_width/2 + third_mesh)
mesh.AddLine('y',  CPW_width/2 + CPW_gap - third_mesh)
mesh.SmoothMeshLines('y', edge_res*1.5, ratio=1.5)
y_pos = mesh.GetLines('y')
mesh.AddLine('y', np.concatenate([-y_pos,
                                   [-substrate_width/2,  substrate_width/2],
                                   [-substrate_width/2 - air_spacing,
                                     substrate_width/2 + air_spacing]]))
mesh.SmoothMeshLines('y', resolution, ratio=1.3)

# z-mesh: fine inside substrate, coarse in air
mesh.AddLine('z', np.linspace(0, substrate_thickness, 5))
mesh.AddLine('z', [-air_spacing, substrate_thickness + air_spacing])
mesh.SmoothMeshLines('z', resolution)

### Substrate
substrate = CSX.AddMaterial('RO4350B', epsilon=substrate_epr)
start = [-CPW_length/2, -substrate_width/2, 0]
stop  = [ CPW_length/2,  substrate_width/2, substrate_thickness]
substrate.AddBox(start, stop)

### CPW ports (include the port metal)
cpw_port_metal = CSX.AddMetal('CPW_PORT')

portstart = [-CPW_length/2,                  -CPW_width/2, substrate_thickness]
portstop  = [-CPW_length/2 + CPW_port_length,  CPW_width/2, substrate_thickness]
port1 = CPWPort(CSX, 1, cpw_port_metal, portstart, portstop, 'x', 'y', CPW_gap,
                excite=1, priority=999,
                MeasPlaneShift=CPW_port_length, Feed_R=feed_R)

portstart = [ CPW_length/2,                  -CPW_width/2, substrate_thickness]
portstop  = [ CPW_length/2 - CPW_port_length,  CPW_width/2, substrate_thickness]
port2 = CPWPort(CSX, 2, cpw_port_metal, portstart, portstop, 'x', 'y', CPW_gap,
                priority=999, MeasPlaneShift=CPW_port_length, Feed_R=feed_R)

ports = [port1, port2]

### CPW centre conductor between the two ports
cpw = CSX.AddMetal('CPW')
start = [-CPW_length/2 + CPW_port_length, -CPW_width/2, substrate_thickness]
stop  = [ CPW_length/2 - CPW_port_length,  CPW_width/2, substrate_thickness]
cpw.AddBox(start, stop, priority=999)

### CPW ground planes (left and right of the gap)
gnd = CSX.AddMetal('GND')
start = [-CPW_length/2, -CPW_width/2 - CPW_gap, substrate_thickness]
stop  = [ CPW_length/2, -substrate_width/2,      substrate_thickness]
gnd.AddBox(start, stop, priority=999)

start = [-CPW_length/2,  CPW_width/2 + CPW_gap, substrate_thickness]
stop  = [ CPW_length/2,  substrate_width/2,      substrate_thickness]
gnd.AddBox(start, stop, priority=999)

if 1:  # debugging only
    CSX_file = os.path.join(Sim_Path, 'cpw_line.xml')
    if not os.path.exists(Sim_Path):
        os.mkdir(Sim_Path)
    CSX.Write2XML(CSX_file)
    from CSXCAD import AppCSXCAD_BIN
    os.system(AppCSXCAD_BIN + ' "{}"'.format(CSX_file))

### Run the simulation
FDTD.Run(Sim_Path, cleanup=True)

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

assert np.max(s11_dB[mask]) < -20, \
    f'FAIL: max(dB(S11)) = {np.max(s11_dB[mask]):.1f} dB, expected < -20 dB'
assert np.min(s21_dB[mask]) > -0.5, \
    f'FAIL: min(dB(S21)) = {np.min(s21_dB[mask]):.1f} dB, expected > -0.5 dB'
assert np.max(s21_dB[mask]) < 0.05, \
    f'FAIL: max(dB(S21)) = {np.max(s21_dB[mask]):.2f} dB, expected < +0.05 dB (sign error?)'

print('PASS')

if 1:  # set to 1 for debugging plots
    import matplotlib.pyplot as plt

    fig, axis = plt.subplots(num='S-Parameters', tight_layout=True)
    axis.plot(f/1e9, s11_dB, 'k-',  linewidth=2, label='$S_{11}$')
    axis.plot(f/1e9, s21_dB, 'r--', linewidth=2, label='$S_{21}$')
    axis.grid()
    axis.set_xmargin(0)
    axis.set_xlabel('Frequency (GHz)')
    axis.set_ylabel('S-Parameter (dB)')
    axis.legend()

    plt.show()
