# -*- coding: utf-8 -*-
"""
 Helical Antenna Tutorial

 Tested with
  - python 3.14
  - openEMS v0.37

 (c) 2015-2025 Thorsten Liebig <thorsten.liebig@gmx.de>
      04-Jan-2026: modified to use matplotlib.pyplot instead of pylab
"""

### Import Libraries
import os, tempfile
import numpy as np

from CSXCAD import CSXCAD

from openEMS import openEMS
from openEMS.results import Results
from openEMS.physical_constants import *

# show the structure before simulating; set True here, or run with SHOW_STRUCTURE=1
show_structure = os.environ.get('SHOW_STRUCTURE', '0') == '1'


### Setup the simulation
Sim_Path = os.path.join(tempfile.gettempdir(), 'Helical_Ant')
force_re_sim = False

unit = 1e-3 # all length in mm

f0 = 2.4e9 # center frequency, frequency of interest!
lambda0 = round(C0/f0/unit) # wavelength in mm
fc = 0.5e9 # 20 dB corner frequency

Helix_radius = 20 # --> diameter is ~ lambda/pi
Helix_turns = 10  # --> expected gain is G ~ 4 * 10 = 40 (16dBi)
Helix_pitch = 30  # --> pitch is ~ lambda/4
Helix_mesh_res = 3

gnd_radius = lambda0/2

# feeding
feed_heigth = 3
feed_R = 120    #feed impedance

# size of the simulation box
SimBox = np.array([1, 1, 1.5])*2.0*lambda0

### Setup FDTD parameter & excitation function
FDTD = openEMS(EndCriteria=1e-4)
FDTD.SetGaussExcite( f0, fc )
FDTD.SetBoundaryCond( ['MUR', 'MUR', 'MUR', 'MUR', 'MUR', 'PML_8'] )

### Setup Geometry & Mesh
CSX = CSXCAD.ContinuousStructure()
FDTD.SetCSX(CSX)
mesh = CSX.GetGrid()
mesh.SetDeltaUnit(unit)

max_res = np.floor(C0 / (f0+fc) / unit / 20) # cell size: lambda/20

# create helix mesh
mesh.AddLine('x', [-Helix_radius, 0, Helix_radius])
mesh.SmoothMeshLines('x', Helix_mesh_res)
# add the air-box
mesh.AddLine('x', [-SimBox[0]/2-gnd_radius,  SimBox[0]/2+gnd_radius])
# create a smooth mesh between specified fixed mesh lines
mesh.SmoothMeshLines('x', max_res, ratio=1.4)

# copy x-mesh to y-direction
mesh.SetLines('y', mesh.GetLines('x'))

# create helix mesh in z-direction
mesh.AddLine('z', [0, feed_heigth, Helix_turns*Helix_pitch+feed_heigth])
mesh.SmoothMeshLines('z', Helix_mesh_res)

# add the air-box
mesh.AddLine('z', [-SimBox[2]/2, max(mesh.GetLines('z'))+SimBox[2]/2 ])
# create a smooth mesh between specified fixed mesh lines
mesh.SmoothMeshLines('z', max_res, ratio=1.4)

### Create the Geometry
## * Create the metal helix using the wire primitive.
## * Create a metal gorund plane as cylinder.
# create a perfect electric conductor (PEC)
helix_metal = CSX.AddMetal('helix' )

ang = np.linspace(0,2*np.pi,21)
coil_x = Helix_radius*np.cos(ang)
coil_y = Helix_radius*np.sin(ang)
coil_z = ang/2/np.pi*Helix_pitch

Helix_x=np.array([])
Helix_y=np.array([])
Helix_z=np.array([])
zpos = feed_heigth
for n in range(Helix_turns-1):
    Helix_x = np.r_[Helix_x, coil_x]
    Helix_y = np.r_[Helix_y, coil_y]
    Helix_z = np.r_[Helix_z ,coil_z+zpos]
    zpos = zpos + Helix_pitch

p = np.array([Helix_x, Helix_y, Helix_z])
helix_metal.AddCurve(p)

# create ground circular ground
gnd = CSX.AddMetal( 'gnd' ) # create a perfect electric conductor (PEC)

# add a box using cylindrical coordinates
start = [0, 0, -0.1]
stop  = [0, 0,  0.1]
gnd.AddCylinder(start, stop, radius=gnd_radius)

# apply the excitation & resist as a current source
start = [Helix_radius, 0, 0]
stop  = [Helix_radius, 0, feed_heigth]
port = FDTD.AddLumpedPort(1 ,feed_R, start, stop, 'z', 1.0, priority=5)

# nf2ff calc
nf2ff = FDTD.CreateNF2FFBox(opt_resolution=[lambda0/15]*3)

### Show the structure (see show_structure above)
if show_structure:
    view = dict(view='iso')
    try:
        from CSXCAD.Viewer import Show
    except ImportError:
        print('Showing the structure needs vtk: pip install vtk')
    else:
        Show(CSX, **view)

### Run the simulation
if force_re_sim or not os.path.exists(os.path.join(Sim_Path, 'et')):
    FDTD.Run(Sim_Path, cleanup=True)

### Postprocessing & plotting
res = Results(FDTD, Sim_Path, f=np.linspace(f0-fc, f0+fc, 501))

## Plot the reflection coefficient S11 above the feed point impedance
fig = res.AddFigure('Feed point', grid=(2, 1), sharex=True,
                    title='reflection coefficient $S_{11}$ and feed point impedance')
fig.AddCurve(res.S(1, 1), axes=0, style='k-', linewidth=2, label='$S_{11}$')
fig.SetAxes(0, ylabel='S11 (dB)')
fig.AddCurve(res.Zin(1), fmt='real', axes=1, style='k-', linewidth=2, label='$\\Re(Z_{in})$')
fig.AddCurve(res.Zin(1), fmt='imag', axes=1, style='r--', linewidth=2, label='$\\Im(Z_{in})$')
fig.SetAxes(1, ylabel='Zin (Ohm)')

## the feed port's voltage above its current, over one time axis
fig = res.AddFigure('Port signals', grid=(2, 1), sharex=True, title='port 1 voltage and current')
fig.AddCurve(res.U(1, domain='t'), axes=0, style='k-', label='$U_1$')
fig.AddCurve(res.I(1, domain='t'), axes=1, style='r-', label='$I_1$')

### Create the NFFF contour
## * the far field at phi=0 degrees, and on a grid over phi for Dmax
theta = np.arange(0., 180., 1.)
phi = np.arange(-180, 180, 2)
print( 'calculating the 3D far field...' )
Dmax_dB = 10*np.log10(res.Dmax(nf2ff.name, f0, theta=theta, phi=phi).GetValue())
cut = res.FarFieldCut(nf2ff.name, f0, phi=0, angles=theta)
D_cut = 10*np.log10(cut.GetXY()[1])
theta_HPBW = theta[ np.where(D_cut<Dmax_dB-3)[0][0] ]

## * Display power and directivity
print('radiated power: Prad = {:.5g} W'.format(res.Prad(nf2ff.name, f0, theta=theta, phi=phi).GetValue()))
print('directivity: Dmax = {:.2f} dBi'.format(Dmax_dB))
print('efficiency: nu_rad = {:.1f} %'.format(100*res.Efficiency(nf2ff.name, f0, 1, theta=theta, phi=phi).GetValue()))
print('theta_HPBW = {:.1f} °'.format(theta_HPBW))

## * Plot the pattern: total and circular polarisations
fig = res.AddFigure('Pattern', title='Frequency: {:.2f} GHz'.format(f0/1e9))
fig.AddCurve(cut, style='k-', linewidth=2, label='|E total|')
fig.AddCurve(res.FarFieldCut(nf2ff.name, f0, phi=0, angles=theta, component='cprh'),
             style='r--', linewidth=2, label='|E CPRH|')
fig.AddCurve(res.FarFieldCut(nf2ff.name, f0, phi=0, angles=theta, component='cplh'),
             style='g-.', linewidth=2, label='|E CPLH|')

# show all plots
res.Show()
