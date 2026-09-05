# -*- coding: utf-8 -*-
"""
Tutorials / 7T MRI Loop Coil

Tested with
  - python 3.14
  - openEMS v0.37

(c) 2013-2026 Thorsten Liebig <thorsten.liebig@gmx.de>

This tutorial models a surface loop coil designed for 7 T MRI (proton Larmor
frequency 298 MHz) placed next to a human head model.  The loop is tuned to
resonance with lumped capacitors.  A pre-converted "Ella" Virtual Family voxel
body model (``Ella_centered_298MHz.h5``) is used when present in the working
directory; otherwise the script falls back to the bundled three-layer
ellipsoidal head phantom (skin / skull / brain, tissue properties at 298 MHz
from the IT'IS database).

"""

### Import Libraries
import os
import warnings
import tempfile
import numpy as np
from matplotlib import pylab as plt

from CSXCAD  import ContinuousStructure
from openEMS import openEMS
from openEMS.physical_constants import *
from openEMS.sar_utils import readSAR
from openEMS.utilities import HDF5Dump


### General Setup
Sim_Path = os.path.join(tempfile.gettempdir(), 'MRI_Loop_Coil')

post_proc_only = False

unit = 1e-3  # all lengths in mm

### Loop Coil Parameters
loop_length  = 80           # length of the loop in z-direction (mm)
loop_width   = 60           # width of the loop in y-direction (mm)
loop_strip_w = 5            # metal strip width (mm)
loop_air_gap = loop_strip_w / 3  # gap width for lumped capacitors (mm)
loop_pos_x   = -130         # x position of the loop plane (mm)
loop_C_gap   = 5.4e-12      # tuning capacitance (F)
loop_port_R  = 2.5          # feed resistance (Ohm)

### Human Body Model Setup
# Pre-converted Ella VF voxel model — create with Convert_VF_DiscMaterial (Octave).
body_model_file   = os.path.join(os.getcwd(), 'Ella_centered_298MHz.h5')
body_model_transform = [
    ('RotateAxis', 'x', np.pi),    # flip head-to-foot
    ('RotateAxis', 'z', np.pi/2),  # rotate nose to +x direction
    ('Translate', [0, 5, -720]),   # centre head at origin (mm)
]

body_box_start = np.array([-120, -150, -200])  # head + shoulder crop (mm)
body_box_stop  = np.array([ 100,  150,  130])

mesh_box_start      = np.array([-120, -80, -120])  # high-res body region (mm)
mesh_box_stop       = np.array([ 100,  80,  120])
mesh_box_resolution = 2   # mm

Air_Box = 200   # air spacer beyond body region (mm)

### FDTD / Excitation Parameters
f0 = 298e6   # center frequency — 7T proton Larmor (Hz)
fc = 300e6   # 20 dB Gaussian corner frequency (Hz)

### Post-Processing Parameters
B1_dyn_range = 3   # decades below the maximum shown in the B1 field maps

### Locate Body Model or Phantom Fallback
use_body_model = os.path.isfile(body_model_file)
if not use_body_model:
    warnings.warn(
        'VF body model not found — using homogeneous ellipsoidal phantom fallback.\n'
        '  Expected: {}'.format(body_model_file)
    )

# Bundled phantom: two levels up from this file → openEMS/resources/phantoms/
_here = os.path.dirname(os.path.abspath(__file__))
phantom_file = os.path.normpath(
    os.path.join(_here, '..', '..', 'resources', 'phantoms', 'phantom_head_298MHz.h5')
)

### FDTD Setup
## * Disabled advanced material cell interpolation and make sure to use an unaveraged constant cell material
## * This is less accurate but is required for SAR averaging according to IEC/IEEE 62704-1
FDTD = openEMS(EndCriteria=1e-4, CellConstantMaterial=True)
FDTD.SetGaussExcite(f0, fc)
FDTD.SetBoundaryCond(['MUR'] * 6)

CSX = ContinuousStructure()
FDTD.SetCSX(CSX)
mesh = CSX.GetGrid()
mesh.SetDeltaUnit(unit)

### Loop Coil Geometry
loop_mat = CSX.AddMetal('loop')
caps_y   = CSX.AddLumpedElement('caps_y', ny='y', C=loop_C_gap)
caps_z   = CSX.AddLumpedElement('caps_z', ny='z', C=loop_C_gap)

x = loop_pos_x         # all loop conductors lie in the plane x = loop_pos_x
W = loop_width  / 2    # half-width  (y)
L = loop_length / 2    # half-length (z)
w = loop_strip_w
g = loop_air_gap

# Horizontal (y-direction) strips at the top and bottom of the loop
loop_mat.AddBox([x, -W,    -L    ], [x, -g/2,  -L + w], priority=10)
loop_mat.AddBox([x, -W,     L - w], [x, -g/2,   L    ], priority=10)
loop_mat.AddBox([x,  g/2,  -L    ], [x,  W,    -L + w], priority=10)
loop_mat.AddBox([x,  g/2,   L - w], [x,  W,     L    ], priority=10)

# Vertical (z-direction) strips on the left and right sides
loop_mat.AddBox([x, -W,    -L + w], [x, -W + w, -g/2 ], priority=10)
loop_mat.AddBox([x, -W,     g/2  ], [x, -W + w,  L - w], priority=10)
loop_mat.AddBox([x,  W - w, -L + w], [x,  W,    -g/2 ], priority=10)
loop_mat.AddBox([x,  W - w,  g/2  ], [x,  W,     L - w], priority=10)

# Three tuning capacitors at the mid-points of the left/right sides and top strip
caps_z.AddBox([x, -W + w/2 - g/2, -g/2], [x, -W + w/2 + g/2, g/2], priority=10)
caps_z.AddBox([x,  W - w/2 - g/2, -g/2], [x,  W - w/2 + g/2, g/2], priority=10)
caps_y.AddBox([x, -g/2, L - w/2 - g/2], [x, g/2, L - w/2 + g/2], priority=10)

# Lumped feed port in the bottom strip gap
port = FDTD.AddLumpedPort(
    port_nr=1, R=loop_port_R,
    start=[loop_pos_x, -g/2, -L + w/2 - g/2],
    stop =[loop_pos_x,  g/2, -L + w/2 + g/2],
    p_dir='y', excite=True
)

### Body Model / Phantom
if use_body_model:
    body_mat = CSX.AddDiscMaterial('body_model', filename=body_model_file, filetype=0, scale=1/unit)
    tr = body_mat.GetTransform()
    for op, *args in body_model_transform:
        tr.AddTransform(op, *args, deg=False)
else:
    body_mat = CSX.AddDiscMaterial('body_model', filename=phantom_file, filetype=0, scale=1/unit)

body_mat.AddBox(body_box_start, body_box_stop, priority=0)

### Mesh Generation
# Seed mesh lines at all loop conductor edges
mesh.AddLine('x', [loop_pos_x])
mesh.AddLine('x', [-mesh_box_resolution / 2, mesh_box_resolution / 2])

for y in [ W, W - w, g/2, -g/2, -(W - w), -W]:
    mesh.AddLine('y', [y])
for z in [ L, L - w, g/2, -g/2, -(L - w), -L]:
    mesh.AddLine('z', [z])

# Body / mesh-box boundaries
mesh.AddLine('x', [mesh_box_start[0], mesh_box_stop[0]])
mesh.AddLine('y', [mesh_box_start[1], mesh_box_stop[1]])
mesh.AddLine('z', [mesh_box_start[2], mesh_box_stop[2]])

# Smooth to mesh_box_resolution inside the body region
mesh.SmoothMeshLines('all', mesh_box_resolution, 1.4)

# Add air spacer and smooth globally (~10 cells / lambda_min)
mesh.AddLine('x', [mesh_box_start[0] - Air_Box, mesh_box_stop[0] + Air_Box])
mesh.AddLine('y', [mesh_box_start[1] - Air_Box, mesh_box_stop[1] + Air_Box])
mesh.AddLine('z', [mesh_box_start[2] - Air_Box, mesh_box_stop[2] + Air_Box])

lambda_min = C0 / (f0 + fc) / unit
mesh.SmoothMeshLines('all', lambda_min / 40, 1.5)

### Field and SAR Dump Boxes
dump_xy_start = body_box_start * np.array([1, 1, 0])
dump_xy_stop  = body_box_stop  * np.array([1, 1, 0])
dump_xz_start = body_box_start * np.array([1, 0, 1])
dump_xz_stop  = body_box_stop  * np.array([1, 0, 1])

hf_xy  = CSX.AddDump('Hf_xy',  dump_type=11, file_type=1, frequency=[f0])
hf_xy.AddBox(dump_xy_start, dump_xy_stop)

sar_xy = CSX.AddDump('SAR_xy', dump_type=20, dump_mode=2, file_type=1, frequency=[f0])
sar_xy.AddBox(dump_xy_start, dump_xy_stop)

hf_xz  = CSX.AddDump('Hf_xz',  dump_type=11, file_type=1, frequency=[f0])
hf_xz.AddBox(dump_xz_start, dump_xz_stop)

sar_xz = CSX.AddDump('SAR_xz', dump_type=20, dump_mode=2, file_type=1, frequency=[f0])
sar_xz.AddBox(dump_xz_start, dump_xz_stop)

### Optional: write XML and launch AppCSXCAD for geometry inspection
if 1:
    if not os.path.exists(Sim_Path):
        os.mkdir(Sim_Path)
    CSX_file = os.path.join(Sim_Path, 'MRI_Loop_Coil.xml')
    CSX.Write2XML(CSX_file)
    from CSXCAD import AppCSXCAD_BIN
    os.system(AppCSXCAD_BIN + ' "{}"'.format(CSX_file))

### Run Simulation
if not post_proc_only:
    FDTD.Run(Sim_Path, cleanup=True)

### Post-Processing

## Port post-processing
freq = np.linspace(f0 - fc, f0 + fc, 501)
port.CalcPort(Sim_Path, freq)

Zin = port.uf_tot / port.if_tot
s11 = port.uf_ref / port.uf_inc
P0_in = float(np.interp(f0, freq, port.P_acc))

## S11 plot
fig, ax = plt.subplots()
ax.plot(freq / 1e6, 20 * np.log10(np.abs(s11)), 'k-', lw=2)
ax.set_xlabel('Frequency (MHz)')
ax.set_ylabel('|S₁₁| (dB)')
ax.set_title('Reflection coefficient S₁₁')
ax.grid(True)

## Admittance plot
fig, ax = plt.subplots()
ax.plot(freq / 1e6, np.real(1. / Zin), 'k-', lw=2, label='real')
ax.plot(freq / 1e6, np.imag(1. / Zin), 'r--', lw=2, label='imag')
ax.set_xlabel('Frequency (MHz)')
ax.set_ylabel('Admittance Y_in (S)')
ax.set_title('Feed port admittance')
ax.legend()
ax.grid(True)

## SAR — axial (xy) and sagittal (xz) planes
# readSAR returns (sar[Nx,Ny,Nz], mesh[m], sar_data); mesh coords are in SI metres
fig, axs = plt.subplots(1, 2, figsize=(12, 5))

sar, sar_mesh, sar_data = readSAR(os.path.join(Sim_Path, 'SAR_xy.h5'))
sar_xy = sar[:, :, 0] / P0_in                          # (Nx, Ny)
X, Y = np.meshgrid(sar_mesh[0] / unit, sar_mesh[1] / unit, indexing='ij')  # m → mm

im_sar_xy = axs[0].pcolormesh(X, Y, sar_xy, shading='auto', cmap='hot')
plt.colorbar(im_sar_xy, ax=axs[0])
axs[0].set_aspect('equal')
axs[0].set_xlabel('x (mm)')
axs[0].set_ylabel('y (mm)')
axs[0].set_title('Local SAR — axial (xy)')

sar, sar_mesh, sar_data = readSAR(os.path.join(Sim_Path, 'SAR_xz.h5'))
sar_xz = sar[:, 0, :] / P0_in                          # (Nx, Nz)
X, Z = np.meshgrid(sar_mesh[0] / unit, sar_mesh[2] / unit, indexing='ij')  # m → mm

im_sar_xz = axs[1].pcolormesh(X, Z, sar_xz, shading='auto', cmap='hot')
plt.colorbar(im_sar_xz, ax=axs[1])
axs[1].set_aspect('equal')
axs[1].set_xlabel('x (mm)')
axs[1].set_ylabel('z (mm)')
axs[1].set_title('Local SAR — sagittal (xz)')

fig.suptitle('SAR / P_in  (W/kg per W)')

# show both SAR plots with the same color range
sar_max = max(sar_xy.max(), sar_xz.max())
im_sar_xy.set_clim(0, sar_max)
im_sar_xz.set_clim(0, sar_max)

## B1 field maps — axial (xy) plane
# the dump is a single xy-plane; SetPlane drops the length-1 z-axis so the
# field comes back as (3, Nx, Ny), with the mesh lines to match, in metres
with HDF5Dump(os.path.join(Sim_Path, 'Hf_xy.h5')) as dump:
    dump.SetPlane('z', pos=0)
    H = dump.GetFieldAtIndex(f_idx=0)
    H_mesh = dump.GetMesh(region=True)
Hx = H[0]
Hy = H[1]
B1p_xy = 0.5 * MUE0 * (Hx + 1j * Hy) / np.sqrt(P0_in)
B1m_xy = 0.5 * MUE0 * (Hx - 1j * Hy) / np.sqrt(P0_in)

X, Y = np.meshgrid(H_mesh['lines'][0] / unit, H_mesh['lines'][1] / unit, indexing='ij')

fig, axs = plt.subplots(1, 2, figsize=(12, 5))
im_B1p_xy = axs[0].pcolormesh(X, Y, np.log10(np.abs(B1p_xy)), shading='auto')
plt.colorbar(im_B1p_xy, ax=axs[0])
axs[0].set_aspect('equal')
axs[0].set_xlabel('x (mm)')
axs[0].set_ylabel('y (mm)')
axs[0].set_title('B₁⁺ field log₁₀ (T/√W) — axial')

im_B1m_xy = axs[1].pcolormesh(X, Y, np.log10(np.abs(B1m_xy)), shading='auto')
plt.colorbar(im_B1m_xy, ax=axs[1])
axs[1].set_aspect('equal')
axs[1].set_xlabel('x (mm)')
axs[1].set_ylabel('y (mm)')
axs[1].set_title('B₁⁻ field log₁₀ (T/√W) — axial')

## B1 field maps — sagittal (xz) plane
with HDF5Dump(os.path.join(Sim_Path, 'Hf_xz.h5')) as dump:
    dump.SetPlane('y', pos=0)
    H = dump.GetFieldAtIndex(f_idx=0)
    H_mesh = dump.GetMesh(region=True)
Hx = H[0]
Hy = H[1]
B1p_xz = 0.5 * MUE0 * (Hx + 1j * Hy) / np.sqrt(P0_in)
B1m_xz = 0.5 * MUE0 * (Hx - 1j * Hy) / np.sqrt(P0_in)

X, Z = np.meshgrid(H_mesh['lines'][0] / unit, H_mesh['lines'][2] / unit, indexing='ij')

fig, axs = plt.subplots(1, 2, figsize=(12, 5))
im_B1p_xz = axs[0].pcolormesh(X, Z, np.log10(np.abs(B1p_xz)), shading='auto')
plt.colorbar(im_B1p_xz, ax=axs[0])
axs[0].set_aspect('equal')
axs[0].set_xlabel('x (mm)')
axs[0].set_ylabel('z (mm)')
axs[0].set_title('B₁⁺ field log₁₀ (T/√W) — sagittal')

im_B1m_xz = axs[1].pcolormesh(X, Z, np.log10(np.abs(B1m_xz)), shading='auto')
plt.colorbar(im_B1m_xz, ax=axs[1])
axs[1].set_aspect('equal')
axs[1].set_xlabel('x (mm)')
axs[1].set_ylabel('z (mm)')
axs[1].set_title('B₁⁻ field log₁₀ (T/√W) — sagittal')

# show all four B1 plots with the same color range,
# covering B1_dyn_range decades below the overall maximum
B1_max = max(np.abs(B1p_xy).max(), np.abs(B1m_xy).max(),
             np.abs(B1p_xz).max(), np.abs(B1m_xz).max())
for im in [im_B1p_xy, im_B1m_xy, im_B1p_xz, im_B1m_xz]:
    im.set_clim(np.log10(B1_max) - B1_dyn_range, np.log10(B1_max))

plt.show()
