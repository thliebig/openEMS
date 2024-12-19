"""
 Microstrip Notch Filter Tutorial
 (c) 2016-2023 Thorsten Liebig <thorsten.liebig@gmx.de>
"""

### Import Libraries
import os, tempfile

# Don't delete this line, we need to test if matplotlib can be imported,
# although it's never used. Incompatible binaries between system and
# PyPI can cause this, see:
# https://github.com/thliebig/openEMS/issues/165
from pylab import *

from CSXCAD  import ContinuousStructure
from openEMS import openEMS
from openEMS.physical_constants import *

### Setup the simulation
Sim_Path = os.path.join(tempfile.gettempdir(), 'NotchFilter')

unit = 1e-6 # specify everything in um
MSL_length = 50000
MSL_width = 600
substrate_thickness = 254
substrate_epr = 3.66
stub_length = 12e3
f_max = 7e9

### Setup FDTD parameters & excitation function

# It's only a smoke test, terminate as soon as possible, otherwise we waste
# CPU time on GitHub Actions (especially for emulated ARM64 systems).
FDTD = openEMS(NrTS=1000)

FDTD.SetGaussExcite( f_max/2, f_max/2 )
FDTD.SetBoundaryCond( ['PML_8', 'PML_8', 'MUR', 'MUR', 'PEC', 'MUR'] )

### Setup Geometry & Mesh
CSX = ContinuousStructure()
FDTD.SetCSX(CSX)
mesh = CSX.GetGrid()
mesh.SetDeltaUnit(unit)

resolution = C0/(f_max*sqrt(substrate_epr))/unit/50 # resolution of lambda/50
third_mesh = array([2*resolution/3, -resolution/3])/4

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

mesh.AddLine('y', [-15*MSL_width, 15*MSL_width+stub_length])
mesh.AddLine('y', (MSL_width/2+stub_length)+third_mesh)
mesh.SmoothMeshLines('y', resolution)

mesh.AddLine('z', linspace(0,substrate_thickness,5))
mesh.AddLine('z', 3000)
mesh.SmoothMeshLines('z', resolution)

## Add the substrate
substrate = CSX.AddMaterial( 'RO4350B', epsilon=substrate_epr)
start = [-MSL_length, -15*MSL_width, 0]
stop  = [+MSL_length, +15*MSL_width+stub_length, substrate_thickness]
substrate.AddBox(start, stop )

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

FDTD.Run(Sim_Path, cleanup=True)

### Post-processing
f = linspace( 1e6, f_max, 1601 )
for p in port:
    p.CalcPort( Sim_Path, f, ref_impedance = 50)

s11 = port[0].uf_ref / port[0].uf_inc
s21 = port[1].uf_ref / port[0].uf_inc
