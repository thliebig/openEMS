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

from CSXCAD  import ContinuousStructure
from openEMS import openEMS
from openEMS.results import Results
from openEMS.physical_constants import *
from openEMS.utilities import HDF5Dump, get_resource_path


### General Setup
Sim_Path = os.path.join(tempfile.gettempdir(), 'MRI_Loop_Coil')

post_proc_only = False
# show the structure before simulating; set True here, or run with SHOW_STRUCTURE=1
show_structure = os.environ.get('SHOW_STRUCTURE', '0') == '1'

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
    # Ellipsoidal head phantom shipped with openEMS
    phantom_file = get_resource_path('phantoms', 'phantom_head_298MHz.h5')

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
        kw = {'deg': False} if op in ('RotateAxis', 'RotateOrigin') else {}
        tr.AddTransform(op, *args, **kw)
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

### Show the structure (see show_structure above)
if show_structure:
    view = dict(view='iso')
    try:
        from CSXCAD.Viewer import Show
    except ImportError:
        print('Showing the structure needs vtk: pip install vtk')
    else:
        Show(CSX, **view)

### Run Simulation
if not post_proc_only:
    FDTD.Run(Sim_Path, cleanup=True)

### Post-Processing

## Port post-processing
res = Results(FDTD, Sim_Path, f=np.linspace(f0 - fc, f0 + fc, 501))
P0_in = float(res.P_acc(1).At(f0))

## S11 above the feed port admittance, over one frequency axis
fig = res.AddFigure('Feed port', grid=(2, 1), sharex=True,
                    title='Reflection coefficient S₁₁ and feed port admittance')
fig.AddCurve(res.S(1, 1), axes=0, style='k-', lw=2, label='$S_{11}$')
fig.SetAxes(0, ylabel='|S₁₁| (dB)')
fig.AddCurve(res.Yin(1), fmt='real', axes=1, style='k-', lw=2, label='real')
fig.AddCurve(res.Yin(1), fmt='imag', axes=1, style='r--', lw=2, label='imag')
fig.SetAxes(1, ylabel='Admittance Y_in (S)')

## the feed port's voltage above its current, over one time axis
fig = res.AddFigure('Port signals', grid=(2, 1), sharex=True, title='port 1 voltage and current')
fig.AddCurve(res.U(1, domain='t'), axes=0, style='k-', label='$U_1$')
fig.AddCurve(res.I(1, domain='t'), axes=1, style='r-', label='$I_1$')

## SAR — axial (xy) and sagittal (xz) planes, per W accepted power, with the same color range
fig = res.AddFigure('SAR', kind='map', grid=(1, 2), shared_clim=True, title='SAR / P_in  (W/kg per W)')
fig.AddMap(res.DumpSlice(os.path.join(Sim_Path, 'SAR_xy.h5'), ('z', 0), freq=f0),
           norm=res.P_acc(1), cmap='hot', axes=0)
fig.AddMap(res.DumpSlice(os.path.join(Sim_Path, 'SAR_xz.h5'), ('y', 0), freq=f0),
           norm=res.P_acc(1), cmap='hot', axes=1)
fig.SetAxes(0, title='Local SAR — axial (xy)', xlabel='x (mm)', ylabel='y (mm)')
fig.SetAxes(1, title='Local SAR — sagittal (xz)', xlabel='x (mm)', ylabel='z (mm)')

## B1 field maps — axial (xy) and sagittal (xz) plane
# the circularly polarised parts of the H field, per sqrt(W) accepted power
def read_B1(fn, plane, pos):
    # the dump is a single plane; SetPlane drops the length-1 axis so the
    # field comes back as (3, N1, N2), with the mesh lines to match, in metres
    with HDF5Dump(os.path.join(Sim_Path, fn)) as dump:
        dump.SetPlane(plane, pos=pos)
        H = dump.GetFieldAtIndex(f_idx=0)
        lines = dump.GetMesh(region=True)['lines']
    B1p = 0.5 * MUE0 * (H[0] + 1j * H[1]) / np.sqrt(P0_in)
    B1m = 0.5 * MUE0 * (H[0] - 1j * H[1]) / np.sqrt(P0_in)
    x, y = [l / unit for n, l in enumerate(lines) if n != 'xyz'.index(plane)]   # m → mm
    return x, y, B1p, B1m

B1 = {'axial': read_B1('Hf_xy.h5', 'z', 0), 'sagittal': read_B1('Hf_xz.h5', 'y', 0)}

# all four B1 plots with the same color range, covering B1_dyn_range decades
# below the overall maximum
B1_max = max(np.abs(b).max() for x, y, *bs in B1.values() for b in bs)
clim = (np.log10(B1_max) - B1_dyn_range, np.log10(B1_max))
for view, (x, y, B1p, B1m) in B1.items():
    fig = res.AddFigure('B1 ' + view, kind='map', grid=(1, 2))
    for axes, (val, name) in enumerate(((B1p, 'B₁⁺'), (B1m, 'B₁⁻'))):
        fig.AddMap(res.ArrayMap(x, y, np.log10(np.abs(val)), label='log₁₀ (T/√W)'), axes=axes, clim=clim)
        fig.SetAxes(axes, title='{} field log₁₀ (T/√W) — {}'.format(name, view), xlabel='x (mm)',
                    ylabel='y (mm)' if view == 'axial' else 'z (mm)')

res.Show()
