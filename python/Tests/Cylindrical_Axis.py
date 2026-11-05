# -*- coding: utf-8 -*-
"""
 Cylindrical coordinates test — E_z on the axis of a closed mesh

 A closed PEC cylinder on a full (closed) cylindrical mesh including the r=0
 axis, excited by two opposite z-directed line sources at r=3 mm on the alpha
 lines pi/4 (+1) and -3pi/4 (-1). The field is odd under a rotation by pi, and
 the mesh is symmetric under it, so E_z on the axis is zero.

 The axis update sums the alpha currents on the ring around the axis; the two
 alpha lines past the seam of a closed mesh repeat its first two lines and
 must not be summed again.

 Pass criteria
   the voltage on the axis stays below 1e-5 of the voltage at a source

 The host extension and the HIP and Metal kernels sum the same currents.
 Engine_Compare.py compares every engine with the basic one on a closed mesh
 including r=0 (its case cylinder_closed).

 Tested with
  - python 3.14
  - openEMS v0.0.37+

"""

import os, tempfile
import numpy as np

from CSXCAD  import ContinuousStructure
from openEMS import openEMS

a_src = (np.pi/4, -3*np.pi/4)

Sim_Path = os.path.join(tempfile.gettempdir(), 'Cylindrical_Axis')

FDTD = openEMS(CoordSystem=1, NrTS=3000, EndCriteria=0)
FDTD.SetGaussExcite(5e9, 4e9)
FDTD.SetBoundaryCond(['PEC'] * 6)   # r=0 and the closed alpha direction need no BC

CSX = ContinuousStructure(CoordSystem=1)
FDTD.SetCSX(CSX)
mesh = CSX.GetGrid()
mesh.SetDeltaUnit(1e-3)
mesh.SetLines('r', np.arange(12))
mesh.SetLines('a', np.linspace(-np.pi, np.pi, 17))
mesh.SetLines('z', np.arange(12) - 5.5)

# a little alpha extent, so the sources snap to their line also where the XML rounds it
e = 1e-4
for name, a, sign in (('p', a_src[0], 1), ('m', a_src[1], -1)):
    CSX.AddExcitation(name, exc_type=0, exc_val=[0, 0, sign]).AddBox([3, a-e, -1.5], [3, a+e, 1.5])
CSX.AddProbe('v_axis', p_type=0).AddBox([0, 0, -0.5], [0, 0, 0.5])
CSX.AddProbe('v_src', p_type=0).AddBox([3, a_src[0], -0.5], [3, a_src[0], 0.5])

FDTD.Run(Sim_Path, cleanup=True)
v_axis = np.max(np.abs(np.loadtxt(os.path.join(Sim_Path, 'v_axis'), comments='%')[:, 1]))
v_src  = np.max(np.abs(np.loadtxt(os.path.join(Sim_Path, 'v_src'), comments='%')[:, 1]))

assert v_src > 0, 'FAIL: the sources produced no field'
print(f'  voltage on the axis {v_axis/v_src:.1e} of the source voltage')
assert v_axis < 1e-5*v_src, \
    f'FAIL: E_z on the axis reaches {v_axis/v_src:.1e} of the source voltage, it must stay zero'

print('PASS')
