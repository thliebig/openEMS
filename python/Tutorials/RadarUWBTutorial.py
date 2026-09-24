# -*- coding: utf-8 -*-
"""
 Tutorials / Radar and UWB: time delay and signal integrity

 Time-domain performance of a simple UWB monopole: the delay from the source
 to the phase centre, and the fidelity -- how similar the radiated waveform
 stays to the excitation. The Gaussian excitation is matched to IEEE 802.15.4
 UWB channel bandwidths, so dispersion and off-resonance effects become
 visible by switching channel.

 Tested with
  - python 3.14
  - openEMS v0.37

 Based on the Octave tutorial by Georg Michel, 2016

 (c) 2026 Thorsten Liebig <thorsten.liebig@gmx.de>

"""

### Import Libraries
import os, tempfile
import numpy as np
import matplotlib.pyplot as plt  # pip install matplotlib

from CSXCAD  import ContinuousStructure
from openEMS import openEMS
from openEMS.physical_constants import *
from openEMS.utilities import DelayFidelity
from openEMS.automesh import mesh_hint_from_box

### Channel Configuration
## Pick one channel. f_c converts the IEEE 802.15.4 20 dB bandwidth to the
## Gaussian 3 dB bandwidth: 3 dB is 0.3925 times the 20 dB bandwidth.
## 'channel4twice' is only there to show the fidelity degrading with bandwidth.
CHANNELS = {
    'channel1'    : (3.5e9, 0.25e9/0.3925),
    'channel2'    : (4.0e9, 0.25e9/0.3925),
    'channel3'    : (4.5e9, 0.25e9/0.3925),
    'channel4'    : (4.0e9, 0.50e9/0.3925),
    'channel5'    : (6.5e9, 0.25e9/0.3925),
    'channel7'    : (6.5e9, 0.50e9/0.3925),
    'channel4twice': (4.0e9, 1.00e9/0.3925),
}
suffix    = 'channel4'
f_0, f_c  = CHANNELS[suffix]

# polarization tilt angle against co-polarization (90 deg is cross polarized)
tilt = 45*np.pi/180

### Antenna and Substrate Parameters
## A planar monopole on FR4 (eps_r = 4). The gap between the feed line and the
## patch sets the impedance match; patchsize is the monopole length.
Sim_Path = os.path.join(tempfile.gettempdir(), 'UWB_' + suffix)
print(f'{Sim_Path=}')

post_proc_only = False
unit = 1e-3     # we will use millimeters

substrate_epsR   = 4      # FR4
substrate_height = 0.707
substrate_cells  = 3      # thickness in cells

gap       = 0.62          # gap between feed line and patch
patchsize = 14

# resolution for the finer structures, e.g. the antenna gap
fineResolution   = C0/(f_0 + f_c)/np.sqrt(substrate_epsR)/unit/40
# resolution for the coarser structures, e.g. the surrounding air
coarseResolution = C0/(f_0 + f_c)/unit/20

### FDTD Configuration
## OverSampling raises the time-domain sample rate well above Nyquist, which
## is what sets the delay resolution of DelayFidelity. All six faces are PML.
FDTD = openEMS(NrTS=30000, EndCriteria=1e-5, OverSampling=20)
FDTD.SetGaussExcite(f_0, f_c)
FDTD.SetBoundaryCond(['PML_8']*6)

CSX = ContinuousStructure()
FDTD.SetCSX(CSX)
mesh = CSX.GetGrid()
mesh.SetDeltaUnit(unit)

### Geometry Setup
## Sheet-like primitives for ground, feed line and patch on a substrate block.
ground    = CSX.AddMetal('Ground')
patch     = CSX.AddMetal('Patch')
line      = CSX.AddMetal('Line')
substrate = CSX.AddMaterial('Substrate', epsilon=substrate_epsR)

substrate.AddBox([-16, -16, -substrate_height], [16, 18, 0], priority=1)
ground.AddBox([-16, -16, -substrate_height], [16, 0, -substrate_height], priority=2)
line_box  = line.AddBox([-1.15, -16, 0], [1.15, gap, 0], priority=2)
patch.AddBox([-patchsize/2, gap, 0], [patchsize/2, gap + patchsize, 0], priority=2)

### Mesh
## Metal edges are added with the 1/3-2/3 rule, then the mesh is graded from
## the fine resolution at the antenna out to the coarse one in free space.
# two mesh lines for the metal coatings of the substrate
mesh.SetLines('z', np.linspace(-substrate_height, 0, substrate_cells + 1))

# edges of the patch and the ground plane, but not yet the microstrip line
FDTD.AddEdges2Grid('xy', properties=patch, metal_edge_res=fineResolution/2)
FDTD.AddEdges2Grid('xy', properties=ground, metal_edge_res=fineResolution/2)

# replace gap mesh lines which are too close by a single mesh line
y = mesh.GetLines('y')
tooclose = np.where(np.diff(y) < fineResolution/4)[0]
if len(tooclose):
    y[tooclose] = (y[tooclose] + y[tooclose+1])/2
    mesh.SetLines('y', np.delete(y, tooclose+1))

# only the x edges of the microstrip, and only the top of the substrate --
# the other substrate edges are covered by the ground plane
FDTD.AddEdges2Grid('x', primitives=line_box, metal_edge_res=fineResolution/2)
mesh.AddLine('y', 18)   # top of the substrate

# so far only edges: now fill in between
mesh.SmoothMeshLines('all', fineResolution)

# add the outer boundary and the coarse free-space lines
mesh.AddLine('x', [-60, 60])
mesh.AddLine('y', [-60, 65])
mesh.AddLine('z', [-46, 45])
mesh.SmoothMeshLines('all', coarseResolution)

### Feeding port and NF2FF box
## The port sits at the outer end of the microstrip, where it reaches the edge
## of the board: the second y edge line of the feed line, i.e. just inside it.
line_hint = mesh_hint_from_box(line_box, 'y', metal_edge_res=fineResolution/2)
feed_y    = sorted(line_hint[1])[1]
y_lines   = mesh.GetLines('y')
feed_y    = y_lines[np.argmin(np.abs(y_lines - feed_y))]   # snap to a mesh line
port = FDTD.AddLumpedPort(1, 50, [-1.15, feed_y, -substrate_height],
                          [1.15, feed_y, 0], 'z', 1.0, priority=999)

x_l, y_l, z_l = (mesh.GetLines(n) for n in 'xyz')
nf2ff = FDTD.CreateNF2FFBox('nf2ff', [x_l[9], y_l[9], z_l[9]],
                            [x_l[-10], y_l[-10], z_l[-10]])

### Run the Simulation
if 0:  # debugging only
    CSX_file = os.path.join(Sim_Path, 'uwb.xml')
    if not os.path.exists(Sim_Path):
        os.makedirs(Sim_Path)
    CSX.Write2XML(CSX_file)
    from CSXCAD import AppCSXCAD_BIN
    os.system(AppCSXCAD_BIN + ' "{}"'.format(CSX_file))

if not post_proc_only:
    FDTD.Run(Sim_Path, cleanup=True)

### Post-Processing
## Reflection coefficient, then delay and fidelity over an azimuthal trace.
freq = np.linspace(f_0 - f_c, f_0 + f_c, 200)
port.CalcPort(Sim_Path, freq)
s11 = port.uf_ref/port.uf_inc
s11phase = np.unwrap(np.angle(s11))

fig, axis = plt.subplots(num="S11", tight_layout=True)
axis.plot(freq/1e6, 20*np.log10(np.abs(s11)), 'k-', linewidth=2)
axis.set_xlabel('frequency f / MHz')
axis.set_ylabel('reflection coefficient $|S_{11}|$ (dB)')
axis.grid()
axis.set_xmargin(0)
axis2 = axis.twinx()
axis2.plot(freq/1e6, s11phase, 'r--', linewidth=2)
axis2.set_ylabel('$S_{11}$ phase (rad)', color='r')
axis2.tick_params(axis='y', colors='r')
axis.set_title(f'reflection coefficient {suffix} $S_{{11}}$')

## Delay and fidelity around the monopole
theta = np.arange(-180, 181, 10)
phi   = [0]
delay, fidelity, res = DelayFidelity(nf2ff, port, Sim_Path,
                                     np.sin(tilt), np.cos(tilt), theta, phi,
                                     f_0, f_c, verbose=1)

# gain at (close to) f_0: turn the directivity into gain
f_idx = int(np.argmin(np.abs(np.asarray(res.freq) - f_0)))
port.CalcPort(Sim_Path, res.freq[f_idx])
gain = res.Dmax[f_idx]*res.Prad[f_idx]/port.P_inc[0]
print(f'gain at {res.freq[f_idx]/1e9:.2f} GHz: {10*np.log10(gain):.2f} dBi')

th = theta/180*np.pi
fig = plt.figure(num="Delay and fidelity", figsize=(11, 5), tight_layout=True)

ax = fig.add_subplot(121, polar=True)
ax.plot(th, delay[:, 0]*C0*1000, 'k-', linewidth=2)
ax.set_theta_zero_location('N')
ax.set_theta_direction(-1)
ax.set_title(f'delay {suffix} / mm')

ax = fig.add_subplot(122, polar=True)
ax.plot(th, fidelity[:, 0]*100, 'r-', linewidth=2)
ax.set_theta_zero_location('N')
ax.set_theta_direction(-1)
ax.set_rlim([98, 100])
ax.set_title(f'fidelity {suffix} / %')

plt.show()
