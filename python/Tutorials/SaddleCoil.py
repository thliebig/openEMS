# -*- coding: utf-8 -*-
"""
Tutorials / NMR Saddle Coil

Tested with
  - python 3.13
  - openEMS v0.37

(c) 2026 B137P107, python port (c) 2026 Sean Mollet

A saddle coil for nuclear magnetic resonance at 500.13 MHz (proton Larmor
frequency at 11.7 T), inside a conducting bore shield. The coil is a single
wire: two axial legs joined by arcs at the top and the bottom, so that the
current in the two halves produces a field transverse to the bore axis.

The run gives the inductance, the resistance and the reflection coefficient
of the coil over the band, and dumps the H-field on two planes for viewing.

A port of the Octave example of https://github.com/thliebig/openEMS/pull/191
"""

### Import Libraries
import os
import tempfile
import numpy as np
from matplotlib import pylab as plt

from CSXCAD import ContinuousStructure
from openEMS import openEMS
from openEMS.physical_constants import *

### General Setup
Sim_Path = os.path.join(tempfile.gettempdir(), 'SaddleCoil')

post_proc_only = False
show_geometry  = False   # open AppCSXCAD instead of simulating

unit = 1e-3   # all lengths in mm

### Excitation
f0 = 500.13e6            # centre frequency (Hz)
fc = 200e6               # 20 dB bandwidth, gaussian excitation

max_res = C0 / (f0 + fc) / unit / 20    # cell size from the highest frequency
Airbox  = C0 / (f0 - fc) / unit / 25    # distance to the boundaries

### Saddle Coil Parameters
coil_height     = 20          # length of the axial legs (mm)
coil_diameter   = 4           # diameter of the coil (mm)
coil_radius     = coil_diameter / 2
coil_angle      = np.pi / 2   # angle the arcs span
wire_radius     = 0.25
wire_space      = 0.25        # gap between the two wire ends at the feed
wire_space_c    = 2 * wire_radius + wire_space
coil_rotation   = coil_angle / 2
coil_copper     = True        # False: perfect electric conductor

port_resist     = 1000        # feed resistance (Ohm)

### Options of the original example
enable_matching_capa = False
enable_tuning_capa   = False
enable_bore_shield   = True

capaM_value = 16e-12
capaT_value = 1.2e-12

### Bore
bore_radius          = 10
bore_shield_thickness = 0.5


def saddle_points():
    """Give the wire path of the coil, one point per row, in mm.

    The path starts at the feed, follows the lower arc, rises along the first
    leg, crosses over at the top and comes back down, so that the current runs
    the same way around both halves.
    """
    ang_space_c = wire_space_c / coil_radius
    pts = []

    def at(angle, radius, z):
        a = angle - coil_rotation
        return [np.cos(a) * radius, np.sin(a) * radius, z]

    # the feed, offset from the coil so the port sits in free space
    a = (np.pi + coil_angle - ang_space_c) / 2 - coil_rotation
    r = coil_radius + wire_space_c
    pts.append([np.cos(a) * r, np.sin(a) * r + wire_space_c, 0])

    for a in np.linspace((np.pi + coil_angle - ang_space_c) / 2,
                         coil_angle - ang_space_c, 9):
        pts.append(at(a, coil_radius + wire_space_c, 0))
    for a in np.linspace(coil_angle - ang_space_c, 0, 9):
        pts.append(at(a, coil_radius, 0))
    for a in np.linspace(0, coil_angle, 9):                      # first leg, up
        pts.append(at(a, coil_radius, coil_height))
    for a in np.linspace(coil_angle, np.pi + coil_angle, 13):    # across the top
        pts.append(at(a, coil_radius, 0))
    for a in np.linspace(np.pi + coil_angle, np.pi, 9):          # second leg
        pts.append(at(a, coil_radius, coil_height))
    pts.append(at(np.pi, coil_radius, wire_space_c))
    pts.append(at(np.pi, coil_radius + wire_space_c, wire_space_c))
    for a in np.linspace(np.pi, (np.pi + coil_angle + ang_space_c) / 2, 9):
        pts.append(at(a, coil_radius + wire_space_c, 0))
    pts.append([pts[-1][0], pts[-1][1] + wire_space_c, pts[-1][2]])

    # cos() and sin() of the same angle do not give the same value twice; round
    # the path so that points which should coincide do
    return np.round(np.array(pts), 3)


### Setup FDTD Parameter & Excitation Function
FDTD = openEMS(NrTS=400000, EndCriteria=1e-5)
FDTD.SetGaussExcite(f0, fc)
FDTD.SetBoundaryCond(['MUR', 'MUR', 'MUR', 'MUR', 'PML_8', 'PML_8'])

CSX = ContinuousStructure()
FDTD.SetCSX(CSX)
mesh = CSX.GetGrid()
mesh.SetDeltaUnit(unit)

### The Coil
points = saddle_points()
if coil_copper:
    coil = CSX.AddMaterial('saddle', kappa=56e6)
else:
    coil = CSX.AddMetal('saddle')
coil.AddWire(points.T, wire_radius, priority=10)

start = points[0]
stop  = points[-1]

if enable_tuning_capa:
    # across the two arc ends, in parallel with the coil
    capaT = CSX.AddLumpedElement('CapaT', ny='x', caps=True, C=capaT_value)
    capaT.AddBox([points[1][0], points[1][1] - 0.25, points[1][2] - 0.25],
                 [points[-2][0], points[-2][1] + 0.25, points[-2][2] + 0.25])

if enable_matching_capa:
    capaM_start = [stop[0] - wire_radius, stop[1], stop[2] - wire_radius]
    capaM_stop  = [stop[0] + wire_radius, stop[1] + wire_space_c, stop[2] + wire_radius]
    capaM = CSX.AddLumpedElement('CapaM', ny='y', caps=True, C=capaM_value)
    capaM.AddBox(capaM_start, capaM_stop)

    # the feed moves to the far side of the matching capacitor
    extend = np.array([start, [start[0], start[1] + wire_space_c, start[2]]])
    coil.AddWire(extend.T, wire_radius, priority=10)
    start = extend[1]
    stop  = [capaM_stop[0] - wire_radius, capaM_stop[1], capaM_stop[2] - wire_radius]

### Excitation Port
port = FDTD.AddLumpedPort(1, port_resist, start, stop, 'x', 1.0, priority=100)

### Independent Voltage and Current Probes
# the port gives u and i as well; these check them against a plain probe
u_probe = CSX.AddProbe('ut1', 0)
u_probe.AddBox(stop, start)

i_probe = CSX.AddProbe('it1', 1)
ref = np.array(start)
i_probe.AddBox([ref[0] - (wire_radius + wire_space / 2), ref[1] - wire_space,
                ref[2] - (wire_radius + wire_space / 2)],
               [ref[0] + (wire_radius + wire_space / 2), ref[1] - wire_space,
                ref[2] + (wire_radius + wire_space / 2)])

### Bore Shield
if enable_bore_shield:
    shield = CSX.AddMetal('BoreShield')
    shield.AddCylindricalShell([0, 0, -wire_radius - Airbox],
                               [0, 0, coil_height + wire_radius + Airbox],
                               bore_radius, bore_shield_thickness, priority=10)

### Mesh
# the edges of the coil and of the shield, then a line every wire radius between
# them, so that the wire keeps its shape in the mesh
edges = {
    'x': [-coil_radius - wire_radius, coil_radius + wire_radius,
          -bore_radius - bore_shield_thickness / 2, bore_radius + bore_shield_thickness / 2],
    'y': [-coil_radius - wire_radius, coil_radius + wire_radius,
          -bore_radius - bore_shield_thickness / 2, bore_radius + bore_shield_thickness / 2],
    'z': [-wire_radius, coil_height + wire_radius],
}
def merge(values, tol):
    """Sort the values and drop those closer together than tol.

    Two arcs of the coil meet at points that differ only in the last digit;
    keeping both would put a cell of a few microns into the mesh and the
    timestep of the whole simulation with it.
    """
    out = []
    for v in sorted(values):
        if not out or v - out[-1] > tol:
            out.append(v)
    return out


for n, ny in enumerate('xyz'):
    lines = merge(np.concatenate([points[:, n], edges[ny]]), wire_radius / 2)
    refined = list(lines)
    for a, b in zip(lines[:-1], lines[1:]):
        N = int(round((b - a) / wire_radius))
        if N >= 2:
            refined.extend(np.linspace(a, b, N))
    refined = merge(refined, wire_radius / 2)
    mesh.AddLine(ny, refined)
    mesh.AddLine(ny, [refined[0] - Airbox, refined[-1] + Airbox])
    mesh.SmoothMeshLines(ny, max_res, 1.4)

# the PML needs its eight cells at the ends of z, at the spacing of the last cell
z = mesh.GetLines('z')
mesh.AddLine('z', [z[0] - (z[1] - z[0]) * n for n in range(1, 9)])
mesh.AddLine('z', [z[-1] + (z[-1] - z[-2]) * n for n in range(1, 9)])

### Field Dumps
# the H-field of the resonance, on the plane across the coil and along the bore
dump_xy = CSX.AddDump('Ht_xy', dump_type=1, file_type=1)
dump_xy.AddBox([-2 * coil_radius, -2 * coil_radius, coil_height / 2],
               [2 * coil_radius, 2 * coil_radius, coil_height / 2])

dump_zy = CSX.AddDump('Ht_zy', dump_type=1, file_type=1)
dump_zy.AddBox([0, -2 * coil_radius, 0], [0, 2 * coil_radius, coil_height])

### Run
if show_geometry:
    CSX_file = os.path.join(Sim_Path, 'saddle.xml')
    if not os.path.exists(Sim_Path):
        os.mkdir(Sim_Path)
    CSX.Write2XML(CSX_file)
    from CSXCAD import AppCSXCAD_BIN
    os.system(AppCSXCAD_BIN + ' "{}"'.format(CSX_file))
elif not post_proc_only:
    FDTD.Run(Sim_Path, verbose=3, cleanup=True)

### Post-Processing
f = np.linspace(f0 - fc, f0 + fc, 501)
port.CalcPort(Sim_Path, f)

Zin = port.uf_tot / port.if_tot
L = np.imag(Zin) / (2 * np.pi * f)
s11 = port.uf_ref / port.uf_inc

idx = np.argmin(np.abs(f - f0))
print('at {:.2f} MHz: inductance {:.1f} nH, resistance {:.2f} Ohm'.format(
      f0 * 1e-6, L[idx] * 1e9, np.real(Zin)[idx]))
# without the tuning and the matching capacitor the coil is an inductor at the
# end of a 1 kOhm feed, so it reflects nearly everything: |S11| says little here
print('|S11| at {:.2f} MHz: {:.1f} dB'.format(f0 * 1e-6, 20 * np.log10(np.abs(s11[idx]))))

plt.figure()
plt.subplot(2, 1, 1)
plt.plot(f * 1e-6, L * 1e9, linewidth=2)
plt.xlabel('Frequency (MHz)')
plt.ylabel('Coil inductance (nH)')
plt.grid()

plt.subplot(2, 1, 2)
plt.plot(f * 1e-6, np.real(Zin), linewidth=2, label='real')
plt.plot(f * 1e-6, np.imag(Zin), 'r', linewidth=2, label='imaginary')
plt.xlabel('Frequency (MHz)')
plt.ylabel('Impedance (Ohm)')
plt.legend()
plt.grid()

plt.figure()
plt.plot(f * 1e-6, 20 * np.log10(np.abs(s11)), 'k-', linewidth=2)
plt.ylim([-40, 10])
plt.title('Reflection coefficient S11')
plt.ylabel('Reflection coefficient |S11| (dB)')
plt.xlabel('Frequency (MHz)')
plt.grid()

plt.show()
