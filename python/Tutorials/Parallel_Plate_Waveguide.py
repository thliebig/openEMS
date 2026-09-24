# -*- coding: utf-8 -*-
"""
 Tutorials / Parallel Plate Waveguide

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

### Setup the simulation
Sim_Path = os.path.join(tempfile.gettempdir(), 'Parallel_Plate_WG')
print(f'{Sim_Path=}')

### FDTD Parameters and Boundary Conditions
## Run 200 time steps with a 10 MHz sinusoidal excitation to reach steady
## state quickly. PEC boundaries on +/-y model the conducting plates; PMC on
## +/-x makes the structure periodic in x; Mur ABCs on +/-z absorb outgoing
## waves.
FDTD = openEMS(NrTS=200, EndCriteria=0, OverSampling=50)
FDTD.SetSinusExcite(10e6)
FDTD.SetBoundaryCond(['PMC', 'PMC', 'PEC', 'PEC', 'MUR', 'MUR'])

### CSXCAD Geometry and Mesh
## All coordinates are in metres. The uniform 1 m mesh spans +/-10 m in x
## and y (the plate aperture) and -10 to 30 m in z, giving 30 cells of
## propagation distance beyond the source plane.
CSX = ContinuousStructure()
FDTD.SetCSX(CSX)
mesh = CSX.GetGrid()
mesh.SetDeltaUnit(1)

mesh.SetLines('x', np.arange(-10, 11, 1))
mesh.SetLines('y', np.arange(-10, 11, 1))
mesh.SetLines('z', np.arange(-10, 31, 1))

### Excitation
## A y-polarised (E_y) uniform-field source at z = 0 launches the TEM
## mode. The excitation box covers the full cross-section to produce a
## spatially uniform plane-wave front.
exc = CSX.AddExcitation('excitation', exc_type=0, exc_val=[0, 1, 0])
exc.AddBox([-10, -10, 0], [10, 10, 0])

### Field Dump
## Record the time-domain E-field in the xz mid-plane (y = 0) so Paraview
## can animate wave propagation along z after the simulation completes.
Et = CSX.AddDump('Et', dump_mode=1)
Et.AddBox([-10, 0, -10], [10, 0, 30])

### Run the simulation
if 0:  # debugging only
    CSX_file = os.path.join(Sim_Path, 'parallel_plate_wg.xml')
    if not os.path.exists(Sim_Path):
        os.mkdir(Sim_Path)
    CSX.Write2XML(CSX_file)
    from CSXCAD import AppCSXCAD_BIN
    os.system(AppCSXCAD_BIN + ' "{}"'.format(CSX_file))

FDTD.Run(Sim_Path, cleanup=True, verbose=3)

print('use Paraview to visualize the FDTD result...')
