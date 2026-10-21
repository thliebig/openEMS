# -*- coding: utf-8 -*-
"""
 Tutorials / Circular Waveguide

 Tested with
  - python 3.14
  - openEMS v0.37

 (c) 2026 Thorsten Liebig <thorsten.liebig@gmx.de>

"""

### Import Libraries
import os, tempfile
import numpy as np

from CSXCAD  import ContinuousStructure
import CSXCAD.Viewer  # early: lets the 3D view open after matplotlib figures
from CSXCAD.CSRectGrid import CoordinateSystem
from openEMS import openEMS
from openEMS.results import Results

### Simulation Parameters
## Define the waveguide dimensions and the frequency range of interest.
## The radius `rad` sets the TE11 cut-off frequency (~251 MHz for 350 mm);
## the chosen frequency span (300-500 MHz) therefore straddles the onset of
## propagation and captures strongly dispersive behaviour.
Sim_Path = os.path.join(tempfile.gettempdir(), 'Circ_WG')
print(f'{Sim_Path=}')

post_proc_only = False
# show the structure before simulating; set True here, or run with SHOW_STRUCTURE=1
show_structure = os.environ.get('SHOW_STRUCTURE', '0') == '1'
# show field results in 3D after simulating; set False here, or run with SHOW_FIELDS=0
show_fields = os.environ.get('SHOW_FIELDS', '1') == '1'
unit = 1e-3    # drawing unit in mm

# waveguide dimensions
length = 2000
rad    = 350   # waveguide radius in mm

# frequency range of interest
f_start = 300e6
f_stop  = 500e6

# targeted mesh resolution (r, a, z)
mesh_res = [10, 2*np.pi/49.999, 10]

### FDTD Parameters and Excitation
## Configure the FDTD solver for cylindrical coordinates (`CoordSystem=1`)
## and drive it with a Gaussian pulse centred in the band. `EndCriteria`
## of 1e-4 stops the time-stepping once the remaining energy falls below
## that fraction of the peak, keeping runtime short while ensuring the
## impulse has decayed completely. PML layers on both z-faces absorb
## outgoing energy without spurious reflections.
FDTD = openEMS(EndCriteria=1e-4, CoordSystem=1)
FDTD.SetGaussExcite(0.5*(f_start+f_stop), 0.5*(f_stop-f_start))

# boundary conditions: pml in pos. and neg. z-direction
FDTD.SetBoundaryCond([0, 0, 0, 0, 3, 3])

### Cylindrical Mesh Setup
## A cylindrical mesh (r, azimuth, z) is the natural coordinate system for
## this geometry: it avoids staircase errors on the curved conducting wall
## and keeps cell counts manageable. SmoothMeshLines distributes lines
## evenly from the axis to the wall in r, over a full 2pi in azimuth, and
## along the waveguide length in z.
CSX = ContinuousStructure(CoordSystem=CoordinateSystem.CYLINDRICAL)
FDTD.SetCSX(CSX)
mesh = CSX.GetGrid()
mesh.SetDeltaUnit(unit)

mesh.SetLines('r', [0, rad])            # mesh in radial direction
mesh.SetLines('a', [0, 2*np.pi])        # mesh in azimuthal direction
mesh.SetLines('z', [0, length])
for ny, res in zip('raz', mesh_res):
    mesh.SmoothMeshLines(ny, res)

r_lines = mesh.GetLines('r')
a_lines = mesh.GetLines('a')
z_lines = mesh.GetLines('z')

### Waveguide Port Definition
## Two TE11 mode ports bookend the waveguide: port 0 (excite=1) injects the
## dominant circular-waveguide mode; port 1 at the far end acts as a passive
## detector. Placing each port several cells inward from the PML boundary
## ensures the mode is fully formed before reaching the absorber and that
## the port plane samples only the travelling wave.
ports = []
start = [r_lines[0],  a_lines[0],  z_lines[7]]
stop  = [r_lines[-1], a_lines[-1], z_lines[14]]
ports.append(FDTD.AddCircWaveGuidePort(0, start, stop, 'z', rad*unit, 'TE11', 0, 1))

start = [r_lines[0],  a_lines[0],  z_lines[-14]]
stop  = [r_lines[-1], a_lines[-1], z_lines[-15]]
ports.append(FDTD.AddCircWaveGuidePort(1, start, stop, 'z', rad*unit, 'TE11'))

### Field Dump Configuration
## Register a volumetric electric-field dump in HDF5 format (`file_type=1`)
## spanning the entire simulation domain. Sub-sampling of 4 in every
## direction reduces the file size by a factor of 64 while still capturing
## the spatial field structure at a resolution sufficient for visualisation
## in ParaView.
Et = CSX.AddDump('Et', file_type=1, sub_sampling=[4, 4, 4])
start = [r_lines[0],  a_lines[0],  z_lines[0]]
stop  = [r_lines[-1], a_lines[-1], z_lines[-1]]
Et.AddBox(start, stop)

### Show the structure (see show_structure above)
if show_structure:
    view = dict(view='iso', grid=dict(outline=True))
    try:
        from CSXCAD.Viewer import Show
    except ImportError:
        print('Showing the structure needs vtk: pip install vtk')
    else:
        Show(CSX, **view)

### Run the simulation
if not post_proc_only:
    FDTD.Run(Sim_Path, cleanup=True)

### Post-Processing: Port Calculation
## CalcPort reads the time-domain probe signals, transforms them to the
## frequency domain, and decomposes incident and reflected wave amplitudes
## at each port. S-parameters follow directly from the voltage wave ratios;
## the complex wave impedance ZL is the ratio of total voltage to total
## current and reveals the dispersive character of the TE11 mode near
## cut-off. Results does this for every quantity of a figure.
res = Results(FDTD, Sim_Path, f=np.linspace(f_start, f_stop, 201))

### S-Parameter Plot
## Plot S11 and S21 in dB across the frequency band. Near cut-off S11
## rises sharply as the mode cannot propagate; well above cut-off S21
## should approach 0 dB (lossless transmission) and S11 should drop,
## confirming that the waveguide is matched to the TE11 mode ports.
fig = res.AddFigure('S-Parameter and ZL', grid=(2, 1), sharex=True)
fig.AddCurve(res.S(0, 0), axes=0, style='k-', linewidth=2, label='$S_{11}$')
fig.AddCurve(res.S(1, 0), axes=0, style='r--', linewidth=2, label='$S_{21}$')

### Waveguide Impedance Comparison
## Overlay the numerically extracted wave impedance (real and imaginary
## parts of ZL) against the analytic TE11 value. Close agreement between
## the two validates both the mode excitation and the port normalisation;
## the strong frequency dependence near cut-off is the hallmark of
## dispersive waveguide propagation.
# below the S-parameters, on the same frequency axis
fig.AddCurve(res.Zin(0), fmt='real', axes=1, style='k-', linewidth=2, label='$\\Re\\{Z_L\\}$')
fig.AddCurve(res.Zin(0), fmt='imag', axes=1, style='r--', linewidth=2, label='$\\Im\\{Z_L\\}$')
fig.AddCurve(res.ZL(0), fmt='real', axes=1, style='g-.', linewidth=2, label='$Z_{L, analytic}$')
fig.SetAxes(1, ylabel='ZL $(\\Omega)$')

### Show the results, with the field in 3D next to the structure (see show_fields above)
## A cut across the guide with the E-field as arrows; play animates the
## recorded time steps.
v = None
if show_fields:
    try:
        from CSXCAD.Viewer import Viewer
        from openEMS.Viewer import FieldSlices
    except ImportError:
        print('Showing the fields needs vtk: pip install vtk')
    else:
        v = Viewer(CSX, view='iso')
        v.AddItem(FieldSlices(os.path.join(Sim_Path, 'Et.h5'), planes=[('z', None)], quantity='magnitude',
                              glyphs=True, glyph_density=15))
if v is None:
    res.Show()
else:
    res.AddToViewer(v)
    v.Show()
