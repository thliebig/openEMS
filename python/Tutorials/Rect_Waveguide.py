# -*- coding: utf-8 -*-
"""
 Rectangular Waveguide Tutorial

 Tested with
  - python 3.14
  - openEMS v0.37

 (c) 2015-2023 Thorsten Liebig <thorsten.liebig@gmx.de>
     04-Jan-2026: modified to use matplotlib.pyplot instead of pylab

"""

### Import Libraries
import os, tempfile
import numpy as np

from CSXCAD  import ContinuousStructure
import CSXCAD.Viewer  # early: lets the 3D view open after matplotlib figures
from openEMS import openEMS
from openEMS.results import Results
from openEMS.physical_constants import *

### Setup the simulation
Sim_Path = os.path.join(tempfile.gettempdir(), 'Rect_WG')

post_proc_only = False
# show the structure before simulating; set True here, or run with SHOW_STRUCTURE=1
show_structure = os.environ.get('SHOW_STRUCTURE', '0') == '1'
# show field results in 3D after simulating; set False here, or run with SHOW_FIELDS=0
show_fields = os.environ.get('SHOW_FIELDS', '1') == '1'
unit = 1e-6; #drawing unit in um

# waveguide dimensions
# WR42
a = 10700;   #waveguide width
b = 4300;    #waveguide height
length = 50000;

# frequency range of interest
f_start = 20e9;
f_0     = 24e9;
f_stop  = 26e9;
lambda0 = C0/f_0/unit;

#waveguide TE-mode definition
TE_mode = 'TE10';

#targeted mesh resolution
mesh_res = lambda0/30

### Setup FDTD parameter & excitation function
FDTD = openEMS(NrTS=1e4);
FDTD.SetGaussExcite(0.5*(f_start+f_stop),0.5*(f_stop-f_start));

# boundary conditions
FDTD.SetBoundaryCond([0, 0, 0, 0, 3, 3]);

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
start=[0, 0, 10*mesh_res];
stop =[a, b, 15*mesh_res];
mesh.AddLine('z', [start[2], stop[2]])
ports.append(FDTD.AddRectWaveGuidePort( 0, start, stop, 'z', a*unit, b*unit, TE_mode, 1))

start=[0, 0, length-10*mesh_res];
stop =[a, b, length-15*mesh_res];
mesh.AddLine('z', [start[2], stop[2]])
ports.append(FDTD.AddRectWaveGuidePort( 1, start, stop, 'z', a*unit, b*unit, TE_mode))

mesh.SmoothMeshLines('all', mesh_res, ratio=1.4)

### Define dump box...
Et = CSX.AddDump('Et', file_type=0, sub_sampling=[2,2,2])
start = [0, 0, 0];
stop  = [a, b, length];
Et.AddBox(start, stop);

if show_fields:
    # the E-field at f_0 as HDF5 frequency-domain dump, for the 3D view below
    Ef = CSX.AddDump('Ef', dump_type=10, file_type=1, dump_mode=2, frequency=[f_0], sub_sampling=[2,2,2])
    Ef.AddBox(start, stop)

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

### Postprocessing & plotting
res = Results(FDTD, Sim_Path, f=np.linspace(f_start, f_stop, 201))

## Plot s-parameter
fig = res.AddFigure('S-Parameter and ZL', grid=(2, 1), sharex=True)
fig.AddCurve(res.S(0, 0), axes=0, style='k-', linewidth=2, label='$S_{11}$')
fig.AddCurve(res.S(1, 0), axes=0, style='r--', linewidth=2, label='$S_{21}$')

## Compare analytic and numerical wave-impedance
# below the S-parameters, on the same frequency axis
fig.AddCurve(res.Zin(0), fmt='real', axes=1, style='k-', linewidth=2, label='$\\Re\\{Z_L\\}$')
fig.AddCurve(res.Zin(0), fmt='imag', axes=1, style='r--', linewidth=2, label='$\\Im\\{Z_L\\}$')
fig.AddCurve(res.ZL(0), fmt='real', axes=1, style='g-.', linewidth=2, label='$Z_{L, analytic}$')
fig.SetAxes(1, ylabel='ZL $(\\Omega)$')

### Show the results, with the field in 3D next to the structure (see show_fields above)
## The slices animate over the phase: press play.
v = None
if show_fields:
    try:
        from CSXCAD.Viewer import Viewer
        from openEMS.Viewer import FieldSlices
    except ImportError:
        print('Showing the fields needs vtk: pip install vtk')
    else:
        v = Viewer(CSX, view='iso_y')
        # the field per volt of the incident port voltage at f_0
        u_inc = res.U(0, part='inc').At(f_0)
        v.AddItem(FieldSlices(os.path.join(Sim_Path, 'Ef.h5'), planes=[('y', b/2), ('x', a/2)],
                              quantity='E_y', norm=u_inc, unit='V/m/V'))
if v is None:
    res.Show()
else:
    res.AddToViewer(v)
    v.Show()
