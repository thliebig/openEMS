# -*- coding: utf-8 -*-
"""
 Tutorials / Conical Horn Antenna

 A conical horn antenna fed by a circular waveguide, excited in its
 dominant TE11 mode.  Horn and feed are built as a single rotationally
 symmetric body from a cross-sectional polygon swept about the z-axis.

 Tested with
  - python 3.14
  - openEMS v0.37

 (c) 2026 Thorsten Liebig <thorsten.liebig@gmx.de>

"""

### Import Libraries
import os, tempfile
import numpy as np

from CSXCAD  import ContinuousStructure
from openEMS import openEMS
from openEMS.results import Results
from openEMS.physical_constants import *

### Simulation Parameters
## Define the length unit and all geometric parameters of the conical horn.
## `horn_radius` is the inner radius of the feeding circular waveguide; its
## value sets the TE11 cut-off frequency and therefore the useful bandwidth
## of the antenna.
Sim_Path = os.path.join(tempfile.gettempdir(), 'Conical_Horn')
print(f'{Sim_Path=}')

post_proc_only = False
# show the structure before simulating; set True here, or run with SHOW_STRUCTURE=1
show_structure = os.environ.get('SHOW_STRUCTURE', '0') == '1'
# show the 3D far-field pattern after simulating; set False here, or run with SHOW_FIELDS=0
show_fields = os.environ.get('SHOW_FIELDS', '1') == '1'
unit = 1e-3     # all length in mm

horn_radius      = 20                # feed waveguide inner radius
horn_length      = 50                # horn length in z-direction
horn_feed_length = 50                # length of the circular feed waveguide
horn_thickness   = 2                 # metal wall thickness
horn_angle       = 20*np.pi/180      # horn opening angle

# size of the simulation box
SimBox = np.array([100, 100, 100])*2

# frequency range of interest
f_start = 10e9
f_stop  = 20e9

# frequency of interest
f0 = 15e9

### FDTD Solver and Excitation Setup
## Initialise the FDTD engine and choose a Gaussian pulse excitation that
## covers the full frequency range of interest in a single simulation run.
## PML absorbing boundaries on all six faces prevent reflections from the
## simulation box edges and emulate an open radiating environment.
FDTD = openEMS(NrTS=30000, EndCriteria=1e-4)
FDTD.SetGaussExcite(0.5*(f_start+f_stop), 0.5*(f_stop-f_start))
FDTD.SetBoundaryCond(['PML_8']*6)

### CSXCAD Geometry and Mesh
## Build the Cartesian mesh that covers both the feeding waveguide (negative
## z) and the radiation half-space above the horn aperture. Fixed lines are
## placed at the waveguide wall and simulation-box boundaries; SmoothMeshLines
## fills in the interior at roughly lambda/15 resolution to resolve the fields
## accurately without oversampling the free-space region.
# currently, openEMS cannot automatically generate a mesh
max_res = C0 / f_stop / unit / 15

CSX = ContinuousStructure()
FDTD.SetCSX(CSX)
mesh = CSX.GetGrid()
mesh.SetDeltaUnit(unit)

# create fixed lines for the simulation box and the waveguide
mesh.SetLines('x', [-SimBox[0]/2, -horn_radius, 0, horn_radius, SimBox[0]/2])
mesh.SmoothMeshLines('x', max_res, ratio=1.4)
mesh.SetLines('y', mesh.GetLines('x'))

mesh.SetLines('z', [-horn_feed_length, 0, SimBox[2]])
mesh.SmoothMeshLines('z', max_res, ratio=1.4)

x_lines = mesh.GetLines('x')
y_lines = mesh.GetLines('y')
z_lines = mesh.GetLines('z')

### Conical Horn Geometry
## Construct the metallic horn and its circular waveguide feed as a single
## rotationally-symmetric body. A cross-sectional polygon is defined in the
## (radial, z) plane and rotated 360 degrees about the z-axis using
## AddRotPoly, which approximates the circular profile on the rectangular
## FDTD grid. The aperture area `A` is precomputed here for use in the
## aperture efficiency calculation during post-processing.
#
# AddRotPoly with norm_dir='x' takes the polygon in the y-z plane, so
# points[0] is the radial distance from the z-axis and points[1] is z.
horn = CSX.AddMetal('Conical_Horn')

r_out = horn_radius + horn_thickness
r_ap  = horn_radius + np.sin(horn_angle)*horn_length   # inner aperture radius
p = np.array([
    [r_out, r_out, r_out + np.sin(horn_angle)*horn_length, r_ap,        horn_radius, horn_radius      ],
    [-horn_feed_length, 0, horn_length,                    horn_length, 0,           -horn_feed_length],
])
horn.AddRotPoly(points=p, norm_dir='x', elevation=0, rot_axis='z',
                angle=[0, 2*np.pi], priority=10)

# horn aperture
A = np.pi * (r_ap*unit)**2

### Waveguide Feed Port
## Excite the dominant TE11 mode inside the circular waveguide using
## AddCircWaveGuidePort. The port spans the lower section of the feed
## waveguide and simultaneously injects the excitation and records the
## incident and reflected wave amplitudes needed to compute S11.
start = [-horn_radius, -horn_radius, z_lines[9]]
stop  = [ horn_radius,  horn_radius, z_lines[0] + horn_feed_length/2]
port = FDTD.AddCircWaveGuidePort(0, start, stop, 'z', horn_radius*unit, 'TE11', 0, 1)

### Excitation Field Dump
## Record a 2-D slice of the electric field at a plane within the feed
## waveguide for visual verification that the TE11 mode profile is being
## launched correctly. Inspecting this dump before analysing results helps
## catch port-placement errors early.
Exc_dump = CSX.AddDump('Exc_dump')
Exc_dump.AddBox([-horn_radius, -horn_radius, z_lines[7]],
                [ horn_radius,  horn_radius, z_lines[7]])

### Near-Field to Far-Field Box
## Place a Huygens surface just inside the PML on five faces to record
## near-field data during the FDTD run. The bottom face (-z direction) is
## excluded because the waveguide feed reaches the lower simulation boundary,
## making that face unsuitable for a Huygens surface.
start = [x_lines[8],  y_lines[8],  z_lines[8]]
stop  = [x_lines[-9], y_lines[-9], z_lines[-9]]
nf2ff = FDTD.CreateNF2FFBox('nf2ff', start, stop, directions=[1, 1, 1, 1, 0, 1])

### Preview Geometry
## Show the structure to inspect the mesh and geometry before committing to
## a full simulation run (see show_structure above). Confirming that the horn
## cross-section and port placement look correct here can save significant
## compute time.
if show_structure:
    view = dict(view='iso', render_mode='solid+edges')
    try:
        from CSXCAD.Viewer import Show
    except ImportError:
        print('Showing the structure needs vtk: pip install vtk')
    else:
        Show(CSX, **view)

### Run Simulation
## Launch the openEMS FDTD solver. The engine iterates time steps until the
## stored energy decays below the EndCriteria threshold, ensuring the
## Fourier-transformed port signals are fully converged before post-processing.
if not post_proc_only:
    FDTD.Run(Sim_Path, cleanup=True)

### Post-Processing and S-Parameter Plot
## Transform the recorded port voltages and currents from the time domain to
## the frequency domain and derive the input reflection coefficient S11.
## A large negative S11 across the band confirms good impedance matching
## between the waveguide feed and the radiating aperture.
res = Results(FDTD, Sim_Path, f=np.linspace(f_start, f_stop, 201))
fig = res.AddFigure('S11', title='reflection coefficient $S_{11}$',
                    ylabel='reflection coefficient $|S_{11}|$ (dB)', ylim=(-60, 0))
fig.AddCurve(res.S(0, 0), style='k-', linewidth=2, label='$S_{11}$')

## the feed port's voltage above its current, over one time axis
fig = res.AddFigure('Port signals', grid=(2, 1), sharex=True, title='port 0 voltage and current')
fig.AddCurve(res.U(0, domain='t'), axes=0, style='k-', label='$U_0$')
fig.AddCurve(res.I(0, domain='t'), axes=1, style='r-', label='$I_0$')

### Far-Field Radiation Patterns
## Invoke the NF2FF transformation at the centre frequency to obtain the
## antenna directivity as a function of elevation angle for two orthogonal
## azimuth cuts (phi = 0 and phi = 90 degrees). The aperture efficiency
## `e_a` relates the achieved directivity to the theoretical maximum for
## a uniformly illuminated aperture of the same physical area.
thetaRange = np.arange(0, 360, 2) - 180
phiRange   = np.array([0, 90])
print('calculating far field at phi=[0 90] deg...')
Dmax = res.Dmax(nf2ff.name, f0, theta=thetaRange, phi=phiRange).GetValue()
Dlog = 10*np.log10(Dmax)
G_a  = 4*np.pi*A/(C0/f0)**2
e_a  = Dmax/G_a

# display some antenna parameter
print(f'radiated power: Prad = {res.Prad(nf2ff.name, f0, theta=thetaRange, phi=phiRange).GetValue()} Watt')
print(f'directivity: Dmax = {Dlog:.4g} dBi')
print(f'aperture efficiency: e_a = {e_a*100:.4g}%')

### Normalised Directivity Plots
## Display the elevation-angle radiation pattern on a linear dB scale and
## as a polar diagram for both azimuth cuts. The polar plot reveals the
## main-lobe beamwidth and side-lobe levels that characterise the antenna's
## angular selectivity and gain performance.
cuts = [(res.FarFieldCut(nf2ff.name, f0, phi=0, angles=thetaRange), 'k-', 'xz-plane ($\\varphi=0^\\circ$)'),
        (res.FarFieldCut(nf2ff.name, f0, phi=90, angles=thetaRange), 'r--', 'yz-plane ($\\varphi=90^\\circ$)')]
fig = res.AddFigure('Pattern', title=f'directivity at {f0/1e9:.0f} GHz', ylim=(-40, 20))
for q, style, label in cuts:
    fig.AddCurve(q, style=style, linewidth=2, label=label)

# polar plot
fig = res.AddFigure('Polar', kind='polar', title=f'directivity (dBi) at {f0/1e9:.0f} GHz',
                    ylim=(-40, 20), theta_zero='N', theta_direction=-1)
for q, style, label in cuts:
    fig.AddCurve(q, style=style, linewidth=2, label=label)

### Three-Dimensional Far-Field Pattern
## Compute the full 3-D radiation pattern by sweeping over a dense grid of
## theta and phi angles. Finer angular spacing is used near the main lobe
## and coarser spacing elsewhere to keep the NF2FF computation fast while
## preserving pattern detail where it matters most. It is shown in 3D next
## to the horn (see show_fields above).
v = None
if show_fields:
    phiRange_3D   = np.unique(np.concatenate((np.arange(-180, -100,  5), np.arange(-100, -50, 2.5),
                                              np.arange( -50,   50,  1), np.arange(  50, 100, 2.5),
                                              np.arange( 100,  181,  5))))
    thetaRange_3D = np.unique(np.concatenate((np.arange(0, 50, 1), np.arange(50, 100, 2),
                                              np.arange(100, 181, 5))))
    print('calculating 3D far field...')
    ff_3D = nf2ff.CalcNF2FF(Sim_Path, f0, thetaRange_3D, phiRange_3D,
                            outfile='nf2ff_3D.h5', verbose=2,
                            read_cached=post_proc_only)
    try:
        from CSXCAD.Viewer import Viewer
        from openEMS.Viewer import FarField
    except ImportError:
        print('Showing the fields needs vtk: pip install vtk')
    else:
        v = Viewer(CSX, view='iso')
        v.AddItem(FarField(ff_3D, quantity='directivity', scale='log'))

# show all plots
if v is None:
    res.Show()
else:
    res.AddToViewer(v)
    v.Show()
