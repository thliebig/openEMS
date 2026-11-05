# -*- coding: utf-8 -*-
"""
 Cylindrical coordinates test — E_z on the axis of a closed mesh

 Two opposite z-directed sources, half a turn apart around the r=0 axis of a
 closed mesh: the field is odd under that rotation, so E_z on the axis must
 stay zero. The axis update sums the alpha currents around it, where the two
 alpha lines past the seam repeat the first two and must not count again.
"""

import os, tempfile
import numpy as np

from CSXCAD  import ContinuousStructure
from openEMS import openEMS

Sim_Path = os.path.join(tempfile.gettempdir(), 'Cylindrical_Axis')

FDTD = openEMS(CoordSystem=1, NrTS=5000, EndCriteria=0)
FDTD.SetGaussExcite(5e9, 4e9)
FDTD.SetBoundaryCond(['PEC'] * 6)
CSX = ContinuousStructure(CoordSystem=1)
FDTD.SetCSX(CSX)
mesh = CSX.GetGrid()
mesh.SetDeltaUnit(1e-3)
mesh.SetLines('r', np.arange(12))
mesh.SetLines('a', np.linspace(-np.pi, np.pi, 17))
mesh.SetLines('z', np.arange(12) - 5.5)
for n, a in enumerate((np.pi/4, -3*np.pi/4)):
    CSX.AddExcitation(f'src{n}', exc_type=0, exc_val=[0, 0, 1-2*n]).AddBox([3, a, -1.5], [3, a, 1.5])
CSX.AddProbe('v_axis', p_type=0).AddBox([0, 0, -0.5], [0, 0, 0.5])
CSX.AddProbe('v_src', p_type=0).AddBox([3, np.pi/4, -0.5], [3, np.pi/4, 0.5])
FDTD.Run(Sim_Path, cleanup=True)

v_axis, v_src = [np.max(np.abs(np.loadtxt(os.path.join(Sim_Path, p), comments='%')[:, 1])) for p in ('v_axis', 'v_src')]
print(f'  voltage on the axis {v_axis/v_src:.1e} of the source voltage')
assert v_axis < 1e-5*v_src, 'FAIL: E_z on the axis must stay zero'
print('PASS')
