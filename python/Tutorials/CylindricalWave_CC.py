# -*- coding: utf-8 -*-
"""
 Tutorials / Cylindrical-Wave Cylindrical Coordinates

 A cylindrical wave launched by an off-centre dipole, simulated on a
 cylindrical mesh with nested azimuthal sub-grids (the "CC" multigrid).

 Tested with
  - python 3.14
  - openEMS v0.37

 (c) 2026 Thorsten Liebig <thorsten.liebig@gmx.de>

"""

### Import Libraries
import os, tempfile
import numpy as np
import matplotlib.pyplot as plt  # pip install matplotlib
from matplotlib.animation import FuncAnimation

from CSXCAD  import ContinuousStructure
from CSXCAD.CSRectGrid import CoordinateSystem
from openEMS import openEMS
from openEMS.utilities import HDF5Dump

### Setup the Simulation
## Define the simulation domain radius, mesh resolution, and five nested
## cylindrical sub-grids whose boundaries progressively double the azimuthal
## cell count, preventing over-sampling of the fields near the axis.
Sim_Path = os.path.join(tempfile.gettempdir(), '2D_CC_Wave')
print(f'{Sim_Path=}')

post_proc_only = False
unit = 1e-3          # drawing unit in mm

mesh_res = 10        # desired mesh resolution
radius   = 2560      # simulation domain radius
split    = [80, 160, 320, 640, 1280]   # radii to split the mesh into sub-grids
split_N  = len(split)                  # number of nested sub-grids
height   = mesh_res*4

f0 = 1e9

excite_offset = 1300
excite_angle  = 45

### FDTD Parameters and Excitation
## CoordSystem=1 selects cylindrical coordinates; MultiGrid activates the
## nested sub-grid engine. A PML on the outer radial face absorbs the
## outgoing cylindrical wave; all other boundaries default to PEC.
FDTD = openEMS(NrTS=100000, EndCriteria=1e-4, CoordSystem=1, MultiGrid=split)
FDTD.SetGaussExcite(f0, f0/2)
FDTD.SetBoundaryCond([0, 3, 0, 0, 0, 0])   # pml in positive r-direction

### CSXCAD Geometry and Mesh
## The outermost sub-domain carries 50 * 2^5 = 1600 azimuthal lines; each
## inner sub-grid halves this count so angular resolution scales with cell
## size. SmoothMeshLines distributes radial and axial lines uniformly.
# 50 mesh lines for the inner most mesh
# increase the total number of meshlines in alpha direction for all sub-grids
N_alpha = 50 * 2**split_N + 1

CSX = ContinuousStructure(CoordSystem=CoordinateSystem.CYLINDRICAL)
FDTD.SetCSX(CSX)
mesh = CSX.GetGrid()
mesh.SetDeltaUnit(unit)

mesh.SetLines('r', [0, radius])
mesh.SmoothMeshLines('r', mesh_res)
mesh.SetLines('a', np.linspace(-np.pi, np.pi, N_alpha))
mesh.SetLines('z', [-height/2, 0, height/2])
mesh.SmoothMeshLines('z', mesh_res)

r_lines = mesh.GetLines('r')
a_lines = mesh.GetLines('a')

### Dipole Excitation
## A z-directed hard E-field excitation placed off-centre at 1300 mm radius
## and 45 degree azimuth launches an asymmetric cylindrical wave, exercising
## the multigrid across its full radial extent.
start = [excite_offset, excite_angle/180*np.pi - 0.001, -20]
stop  = [excite_offset, excite_angle/180*np.pi + 0.001,  20]
if excite_offset == 0:
    start[1] = a_lines[0]
    stop[1]  = a_lines[0]

exc = CSX.AddExcitation('excite', exc_type=1, exc_val=[0, 0, 1])
exc.AddBox(start, stop)

### Field Dump Boxes
## Two overlapping dump regions cover the full r-alpha plane at z = 0.
## The time-domain VTK dump is sub-sampled for Paraview; the frequency-domain
## HDF5 dump stores the complex E-field phasor at f0 for post-processing.
start = [r_lines[0],  a_lines[0],  0]
stop  = [r_lines[-9], a_lines[-1], 0]

# time domain vtk dump
Et = CSX.AddDump('Et_ra', dump_type=0, file_type=0, sub_sampling=[4, 10, 1])
Et.AddBox(start, stop)

# frequency domain hdf5 dump
Ef = CSX.AddDump('Ef_ra', dump_type=10, file_type=1, sub_sampling=[2, 2, 2])
Ef.SetFrequency([f0])
Ef.AddBox(start, stop)

### Run the simulation
if 0:  # debugging only
    CSX_file = os.path.join(Sim_Path, '2D_CC_Wave.xml')
    if not os.path.exists(Sim_Path):
        os.mkdir(Sim_Path)
    CSX.Write2XML(CSX_file)
    from CSXCAD import AppCSXCAD_BIN
    os.system(AppCSXCAD_BIN + ' "{}"'.format(CSX_file))

if not post_proc_only:
    FDTD.Run(Sim_Path, cleanup=True)

### Paraview Visualization
## The time-domain VTK dump can be opened in Paraview to animate the
## propagating wave front directly on the cylindrical mesh.
print('use Paraview to visualize the vtk field dump...')

### Post-processing and Phase Animation
## Read the frequency-domain HDF5 dump, convert the cylindrical mesh to
## Cartesian coordinates, then animate the E_z phasor over 0-360 degrees
## to visualise the full cylindrical wave pattern.
with HDF5Dump(os.path.join(Sim_Path, 'Ef_ra.h5')) as dump:
    h5_mesh = dump.GetMesh()
    Ez = dump.GetFieldAtFrequency(f0, component='z')

r = h5_mesh['lines'][0]
a = h5_mesh['lines'][1]

a  = np.append(a, a[0])             # closeup mesh for visualization
Ez = np.squeeze(Ez)
Ez = np.concatenate((Ez, Ez[:, :1]), axis=1)

R, A = np.meshgrid(r, a, indexing='ij')
X = R*np.cos(A)
Y = R*np.sin(A)

E_max = np.max(np.abs(Ez))          # get maximum E_z amplitude

fig, axis = plt.subplots(num="Ez", tight_layout=True)
quad = axis.pcolormesh(X, Y, np.real(Ez), cmap='RdBu_r',
                       vmin=-E_max/10, vmax=E_max/10, shading='gouraud')
axis.set_aspect('equal')
axis.set_xlabel('x (m)')
axis.set_ylabel('y (m)')
axis.set_title(f'$E_z$ at {f0/1e9:.1f} GHz')
fig.colorbar(quad, ax=axis, label='$E_z$ (V/m)')

def _phase(ph):
    """animate phase from 0..360 degree"""
    quad.set_array(np.real(Ez*np.exp(1j*ph/180*np.pi)))
    return (quad,)

# keep a reference, otherwise the animation is garbage collected
anim = FuncAnimation(fig, _phase, frames=np.linspace(0, 360, 41),
                     interval=100, blit=False)

plt.show()
