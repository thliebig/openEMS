# -*- coding: utf-8 -*-
"""
 Simple Patch Antenna Tutorial

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

### General parameter setup
Sim_Path = os.path.join(tempfile.gettempdir(), 'Simp_Patch')

post_proc_only = False
# show the structure before simulating; set True here, or run with SHOW_STRUCTURE=1
show_structure = os.environ.get('SHOW_STRUCTURE', '0') == '1'
# show field results in 3D after simulating; set False here, or run with SHOW_FIELDS=0
show_fields = os.environ.get('SHOW_FIELDS', '1') == '1'

# patch width (resonant length) in x-direction
patch_width  = 32 #
# patch length in y-direction
patch_length = 40

#substrate setup
substrate_epsR   = 3.38
substrate_kappa  = 1e-3 * 2*np.pi*2.45e9 * EPS0*substrate_epsR
substrate_width  = 60
substrate_length = 60
substrate_thickness = 1.524
substrate_cells = 4

#setup feeding
feed_pos = -6 #feeding position in x-direction
feed_R = 50     #feed resistance

# size of the simulation box
SimBox = np.array([200, 200, 150])

# setup FDTD parameter & excitation function
f0 = 2e9 # center frequency
fc = 1e9 # 20 dB corner frequency

### FDTD setup
## * Limit the simulation to 30k timesteps
## * Define a reduced end criteria of -40dB
FDTD = openEMS(NrTS=30000, EndCriteria=1e-4)
FDTD.SetGaussExcite( f0, fc )
FDTD.SetBoundaryCond( ['MUR', 'MUR', 'MUR', 'MUR', 'MUR', 'MUR'] )


CSX = ContinuousStructure()
FDTD.SetCSX(CSX)
mesh = CSX.GetGrid()
mesh.SetDeltaUnit(1e-3)
mesh_res = C0/(f0+fc)/1e-3/20

### Generate properties, primitives and mesh-grid
#initialize the mesh with the "air-box" dimensions
mesh.AddLine('x', [-SimBox[0]/2, SimBox[0]/2])
mesh.AddLine('y', [-SimBox[1]/2, SimBox[1]/2]          )
mesh.AddLine('z', [-SimBox[2]/3, SimBox[2]*2/3]        )

# create patch
patch = CSX.AddMetal( 'patch' ) # create a perfect electric conductor (PEC)
start = [-patch_width/2, -patch_length/2, substrate_thickness]
stop  = [ patch_width/2 , patch_length/2, substrate_thickness]
patch.AddBox(priority=10, start=start, stop=stop) # add a box-primitive to the metal property 'patch'
FDTD.AddEdges2Grid(dirs='xy', properties=patch, metal_edge_res=mesh_res/2)

# create substrate
substrate = CSX.AddMaterial( 'substrate', epsilon=substrate_epsR, kappa=substrate_kappa)
start = [-substrate_width/2, -substrate_length/2, 0]
stop  = [ substrate_width/2,  substrate_length/2, substrate_thickness]
substrate.AddBox( priority=0, start=start, stop=stop )

# add extra cells to discretize the substrate thickness
mesh.AddLine('z', np.linspace(0,substrate_thickness,substrate_cells+1))

# create ground (same size as substrate)
gnd = CSX.AddMetal( 'gnd' ) # create a perfect electric conductor (PEC)
start[2]=0
stop[2] =0
gnd.AddBox(start, stop, priority=10)

FDTD.AddEdges2Grid(dirs='xy', properties=gnd)

# apply the excitation & resist as a current source
start = [feed_pos, 0, 0]
stop  = [feed_pos, 0, substrate_thickness]
port = FDTD.AddLumpedPort(1, feed_R, start, stop, 'z', 1.0, priority=5, edges2grid='xy')

mesh.SmoothMeshLines('all', mesh_res, 1.4)

# Add the nf2ff recording box
nf2ff = FDTD.CreateNF2FFBox()

if show_fields:
    # the E-field in time domain as HDF5 dump, in a plane inside the substrate
    # just below the patch; shown at the resonance found below
    # cell interpolation: at the top surface it gives the centre of the top
    # substrate cell, with all field components at the same point
    Et = CSX.AddDump('Et', file_type=1, dump_mode=2)
    z_dump = substrate_thickness
    Et.AddBox([-substrate_width/2, -substrate_length/2, z_dump], [substrate_width/2, substrate_length/2, z_dump])

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
# what to compute and plot; frequencies: the excitation band f0 +- fc
res = Results(FDTD, Sim_Path)

# the input matching above the input impedance, over one frequency axis
fig = res.AddFigure('Feed point', grid=(2, 1), sharex=True, title='Input matching and impedance')
fig.AddCurve(res.S(1, 1), fmt='dB', axes=0, style='k-', linewidth=2, label='$S_{11}$')
fig.AddCurve(res.Zin(1), fmt='real', axes=1, style='k-', linewidth=2, label='$\\Re\\{Z_{in}\\}$')
fig.AddCurve(res.Zin(1), fmt='imag', axes=1, style='r--', linewidth=2, label='$\\Im\\{Z_{in}\\}$')

## the feed port's voltage above its current, over one time axis
fig = res.AddFigure('Port signals', grid=(2, 1), sharex=True, title='port 1 voltage and current')
fig.AddCurve(res.U(1, domain='t'), axes=0, style='k-', label='$U_1$')
fig.AddCurve(res.I(1, domain='t'), axes=1, style='r-', label='$I_1$')

field_items = []   # 3D field displays, for the viewer below

# the numbers behind a curve
f, s11 = res.S(1, 1).GetXY()
s11_dB = 20.0*np.log10(np.abs(s11))

idx = np.where((s11_dB<-10) & (s11_dB==np.min(s11_dB)))[0]
if not len(idx)==1:
    print('No resonance frequency found for far-field calulation')
else:
    f_res = f[idx[0]]
    theta = np.arange(-180.0, 180.0, 2.0)
    fig = res.AddFigure('Pattern', title='Frequency: {} GHz'.format(f_res/1e9), xlabel='Theta (deg)',
                        ylabel='Directivity (dBi)')
    fig.AddCurve(res.FarFieldCut(nf2ff.name, f_res, phi=0, angles=theta, center=[0,0,1e-3]),
                 style='k-', linewidth=2, label='xz-plane')
    fig.AddCurve(res.FarFieldCut(nf2ff.name, f_res, phi=90, angles=theta, center=[0,0,1e-3]),
                 style='r--', linewidth=2, label='yz-plane')

    ### The 3D far-field pattern and the field below the patch (see show_fields above)
    if show_fields:
        # a 3D pattern needs a full theta/phi grid
        nf2ff_3D = nf2ff.CalcNF2FF(Sim_Path, f_res, np.arange(0, 181, 3), np.arange(0, 361, 5),
                                   center=[0,0,1e-3], outfile='nf2ff_3D.h5')
        try:
            from openEMS.Viewer import FarField, FieldSlices
        except ImportError:
            pass   # no vtk: the figures are shown with matplotlib below
        else:
            # E_z below the patch at the resonance (a DFT of the time-domain
            # dump) per volt of the incident port voltage, animated over the
            # phase: press play
            u_inc = res.U(1, part='inc').At(f_res)
            field_items.append(FieldSlices(os.path.join(Sim_Path, 'Et.h5'), quantity='E_z', domain='frequency',
                                           freq=f_res, warp=True, contours=True, legend=False,
                                           norm=u_inc, unit='V/m/V'))
            # the pattern is moved up along its main lobe to clear the field
            # (origin='auto')
            field_items.append(FarField(nf2ff_3D, quantity='directivity', scale='log'))


# store the setup and the results next to the simulation data, to show them
# again later without this script:
#   FDTD = openEMS(); FDTD.ReadFromXML(xml); FDTD.ReadSetup(json)
#   Results.Load(FDTD, Sim_Path, json).Show()
FDTD.Write2XML(os.path.join(Sim_Path, 'Simple_Patch_Antenna.xml'))
res.Write(os.path.join(Sim_Path, 'Simple_Patch_Antenna.json'))

### Show the results, with the fields in 3D next to the structure (see
# show_fields above)
if not field_items:
    res.Show()
else:
    from CSXCAD.Viewer import Viewer
    # only the see-through patch, no mesh grid
    view = dict(
        size=(0.489, 1.047),
        transparency=0.5,
        grid={'visible': False},
        isolate=('patch',),
        camera={'position': (116.541, 163.599, 131.925), 'focal_point': (-4.56237, 5.71686, 35.0712), 'view_up': (-0.302527, -0.318947, 0.898193), 'view_angle': 30.0, 'parallel_scale': 69.3971, 'parallel_projection': False},
    )
    v = Viewer(CSX, **view)
    for item in field_items:
        v.AddItem(item)
    res.AddToViewer(v)
    v.Show()
