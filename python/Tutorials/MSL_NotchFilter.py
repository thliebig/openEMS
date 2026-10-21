# -*- coding: utf-8 -*-
"""
 Microstrip Notch Filter Tutorial

 Tested with
  - python 3.14
  - openEMS v0.37

 (c) 2016-2023 Thorsten Liebig <thorsten.liebig@gmx.de>
 15-Dec-2025: modified to use matplotlib.pyplot instead of pylab

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
Sim_Path = os.path.join(tempfile.gettempdir(), 'NotchFilter')
post_proc_only = False
# show the structure before simulating; set True here, or run with SHOW_STRUCTURE=1
show_structure = os.environ.get('SHOW_STRUCTURE', '0') == '1'
# show field results in 3D after simulating; set False here, or run with SHOW_FIELDS=0
show_fields = os.environ.get('SHOW_FIELDS', '1') == '1'

unit = 1e-6 # specify everything in um
MSL_length = 50000
MSL_width = 600
substrate_thickness = 254
substrate_epr = 3.66
stub_length = 12e3
f_max = 7e9

### Setup FDTD parameters & excitation function
FDTD = openEMS()
FDTD.SetGaussExcite( f_max/2, f_max/2 )
FDTD.SetBoundaryCond( ['PML_8', 'PML_8', 'MUR', 'MUR', 'PEC', 'MUR'] )

### Setup Geometry & Mesh
CSX = ContinuousStructure()
FDTD.SetCSX(CSX)
mesh = CSX.GetGrid()
mesh.SetDeltaUnit(unit)

resolution = C0/(f_max*np.sqrt(substrate_epr))/unit/50 # resolution of lambda/50
third_mesh = np.array([2*resolution/3, -resolution/3])/4

## Do manual meshing
mesh.AddLine('x', 0)
mesh.AddLine('x',  MSL_width/2+third_mesh)
mesh.AddLine('x', -MSL_width/2-third_mesh)
mesh.SmoothMeshLines('x', resolution/4)

mesh.AddLine('x', [-MSL_length, MSL_length])
mesh.SmoothMeshLines('x', resolution)

mesh.AddLine('y', 0)
mesh.AddLine('y',  MSL_width/2+third_mesh)
mesh.AddLine('y', -MSL_width/2-third_mesh)
mesh.SmoothMeshLines('y', resolution/4)

mesh.AddLine('y', (MSL_width/2+stub_length)+third_mesh)
mesh.AddLine('y', [-15*MSL_width, 15*MSL_width+stub_length])
mesh.SmoothMeshLines('y', resolution)

mesh.AddLine('z', np.linspace(0,substrate_thickness,5))
mesh.AddLine('z', 3000)
mesh.SmoothMeshLines('z', resolution)

## Add the substrate
substrate = CSX.AddMaterial( 'RO4350B', epsilon=substrate_epr)
start = [mesh.GetLine('x', 0) , mesh.GetLine('y', 0) , 0]
stop  = [mesh.GetLine('x', -1), mesh.GetLine('y', -1), substrate_thickness]
substrate.AddBox(start, stop)

## MSL port setup
port = [None, None]
pec = CSX.AddMetal( 'PEC' )
portstart = [ -MSL_length, -MSL_width/2, substrate_thickness]
portstop  = [ 0,  MSL_width/2, 0]
port[0] = FDTD.AddMSLPort( 1,  pec, portstart, portstop, 'x', 'z', excite=-1, FeedShift=10*resolution, MeasPlaneShift=MSL_length/3, priority=10)

portstart = [MSL_length, -MSL_width/2, substrate_thickness]
portstop  = [0         ,  MSL_width/2, 0]
port[1] = FDTD.AddMSLPort( 2, pec, portstart, portstop, 'x', 'z', MeasPlaneShift=MSL_length/3, priority=10 )

## Filter-Stub Definition
start = [-MSL_width/2,  MSL_width/2, substrate_thickness]
stop  = [ MSL_width/2,  MSL_width/2+stub_length, substrate_thickness]
pec.AddBox(start, stop, priority=10 )

if show_fields:
    # the E-field in time domain as HDF5 dump, in a plane inside the substrate
    # just below the traces, for the 3D view below
    # cell interpolation: at the top surface it gives the centre of the top
    # substrate cell, with all field components at the same point
    Et = CSX.AddDump('Et', file_type=1, dump_mode=2)
    z_dump = substrate_thickness
    Et.AddBox([mesh.GetLine('x', 0), mesh.GetLine('y', 0), z_dump], [mesh.GetLine('x', -1), mesh.GetLine('y', -1), z_dump])

### Show the structure (see show_structure above)
if show_structure:
    view = dict(view='iso', render_mode='solid+edges')
    try:
        from CSXCAD.Viewer import Show
    except ImportError:
        print('Showing the structure needs vtk: pip install vtk')
    else:
        Show(CSX, **view)

### Run the simulation
if not post_proc_only:
    FDTD.Run(Sim_Path, cleanup=True)

### Post-processing and plotting
res = Results(FDTD, Sim_Path, f=np.linspace(1e6, f_max, 1601))
for p in port:
    res.SetPort(p.number, ref_impedance=50)

fig = res.AddFigure('S-Parameter')
fig.AddCurve(res.S(1, 1), style='k-', linewidth=2, label='$S_{11}$')
fig.AddCurve(res.S(2, 1), style='r--', linewidth=2, label='$S_{21}$')

### Show the results, with the field in 3D next to the structure (see show_fields above)
## The pulse travels along the line and into the stub: press play.
v = None
if show_fields:
    try:
        from CSXCAD.Viewer import Viewer
        from openEMS.Viewer import FieldSlices
    except ImportError:
        print('Showing the fields needs vtk: pip install vtk')
    else:
        # only the see-through traces (to see the field below them), no mesh grid
        v = Viewer(CSX, view='iso', transparency=0.5, isolate='PEC', grid=dict(visible=False))
        v.AddItem(FieldSlices(os.path.join(Sim_Path, 'Et.h5'), quantity='E_z', warp=True))
if v is None:
    res.Show()
else:
    res.AddToViewer(v)
    v.Show()
