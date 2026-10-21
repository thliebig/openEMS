# -*- coding: utf-8 -*-
"""
 Tutorials / radar cross section of a metal sphere

 Tested with
  - python 3.14
  - openEMS v0.37

 (c) 2016-2023 Thorsten Liebig <thorsten.liebig@gmx.de>
"""

### Import Libraries
import os, tempfile
from pylab import *

from CSXCAD  import ContinuousStructure
from openEMS import openEMS
from openEMS.physical_constants import *
from openEMS.results import Results

### Setup the simulation
Sim_Path = os.path.join(tempfile.gettempdir(), 'RCS_Sphere')
post_proc_only = False
# show the structure before simulating; set True here, or run with SHOW_STRUCTURE=1
show_structure = os.environ.get('SHOW_STRUCTURE', '0') == '1'

unit = 1e-3 # all length in mm

sphere_rad = 200

inc_angle = 0 #incident angle (to x-axis) in deg

# size of the simulation box
SimBox = 1200
PW_Box = 750

### Setup FDTD parameters & excitation function
FDTD = openEMS(EndCriteria=1e-5)

f_start =  50e6 # start frequency
f_stop  = 1000e6 # stop  frequency
f0      = 500e6
FDTD.SetGaussExcite(  0.5*(f_start+f_stop), 0.5*(f_stop-f_start) )

FDTD.SetBoundaryCond( ['PML_8', 'PML_8', 'PML_8', 'PML_8', 'PML_8', 'PML_8'] )

### Setup Geometry & Mesh
CSX = ContinuousStructure()
FDTD.SetCSX(CSX)
mesh = CSX.GetGrid()
mesh.SetDeltaUnit(unit)

#create mesh
mesh.SetLines('x', [-SimBox/2, 0, SimBox/2])
mesh.SmoothMeshLines('x', C0 / f_stop / unit / 20) # cell size: lambda/20
mesh.SetLines('y', mesh.GetLines('x'))
mesh.SetLines('z', mesh.GetLines('x'))

### Create a metal sphere and plane wave source
sphere_metal = CSX.AddMetal( 'sphere' ) # create a perfect electric conductor (PEC)
sphere_metal.AddSphere(priority=10, center=[0, 0, 0], radius=sphere_rad)

# plane wave excitation
k_dir = [cos(np.deg2rad(inc_angle)), sin(np.deg2rad(inc_angle)), 0] # plane wave direction
E_dir = [0, 0, 1] # plane wave polarization --> E_z

pw_exc = CSX.AddExcitation('plane_wave', exc_type=10, exc_val=E_dir)
pw_exc.SetPropagationDir(k_dir)
pw_exc.SetFrequency(f0)

start = np.array([-PW_Box/2, -PW_Box/2, -PW_Box/2])
stop  = -start
pw_exc.AddBox(start, stop)

# nf2ff calc
nf2ff = FDTD.CreateNF2FFBox()

### Show the structure (see show_structure above)
if show_structure:
    view = dict(view='iso')
    try:
        from CSXCAD.Viewer import Show
    except ImportError:
        print('Showing the structure needs vtk: pip install vtk')
    else:
        pw_exc.SetRenderMode('edges')  # see the sphere inside the excitation box
        Show(CSX, **view)

### Run the simulation
if not post_proc_only:
    FDTD.Run(Sim_Path, cleanup=True)

### Postprocessing & plotting
# the incident field is the excitation signal 'et' of the plane wave
res = Results(FDTD, Sim_Path)
E_amp = norm(E_dir)

## the bistatic RCS at f0 in the xy-plane
fig = res.AddFigure('RCS pattern', kind='polar', title='RCS at {:g} MHz'.format(f0/1e6))
fig.AddCurve(res.RCS(nf2ff.name, f0, 'et', theta=90, angles=arange(-180, 180.1, 2), E_amp=E_amp),
             style='k-', linewidth=2, label='xy-plane')

## the backscatter over frequency
freq = linspace(f_start, f_stop, 100)
back = res.RCS(nf2ff.name, freq, 'et', theta=90, phi=180+inc_angle, E_amp=E_amp)
fig = res.AddFigure('RCS', title='radar cross section')
fig.AddCurve(back, linewidth=2, label='backscatter')

## normalized to the sphere's cross section, over its radius in wavelengths
f, back_scat = back.GetXY()
a = sphere_rad*unit
fig = res.AddFigure('normalized RCS', title='normalized radar cross section',
                    xlabel='sphere radius / wavelength', ylabel='RCS / ($\\pi a^2$)',
                    yscale='log', ylim=(0.01, 10))
fig.AddCurve(res.Array(a/C0*f, back_scat/(pi*a**2), domain=None), linewidth=2, label='backscatter')

res.Show()
