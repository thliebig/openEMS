# -*- coding: utf-8 -*-
"""
 Sinusoidal excitation period test

 The operator shrinks the timestep of a sinusoidal excitation to a whole
 number of timesteps per period, and the excitation repeats its signal table
 with that period.  Here the period divided by the shrunk timestep comes out
 just below that whole number (2.4 GHz on this mesh: 433 timesteps), the
 case a truncated period gets wrong.  A TEM plane wave (Ex, travelling in +z)
 in a 1D channel (PEC walls normal to x, PMC walls normal to y, PML at both
 z-ends) is driven until the steady-state detection ends the run.

 Pass criteria
   openEMS runs the period in 433 timesteps
   the run ends by steady-state detection (end time < 50 excitation periods)
   the probe signal repeats with the excitation period: its phasor changes
     < 0.1 % between the last two periods

 Tested with
  - python 3.13
  - openEMS v0.0.37+

"""

import os, re, tempfile
import numpy as np

from CSXCAD  import ContinuousStructure
from openEMS import openEMS

### Geometry
unit     = 1e-3    # drawing unit: mm
dz       = 0.5     # axial mesh resolution (mm)
width    = 1.0     # transverse channel width (mm)
length   = 150     # channel length (mm)
z_src    = 20      # excitation plane
z_probe  = 110     # field probe

### Excitation
f0      = 2.4e9
T0      = 1 / f0
TS0     = 433      # timesteps per period on this mesh
NrTS    = 200*TS0  # the steady-state detection must end the run long before

# the case of this test: truncated, the period would come out one timestep short
assert int(T0/(T0/TS0)) == TS0 - 1

Sim_Path = os.path.join(tempfile.gettempdir(), 'Sinusoidal_Period')

FDTD = openEMS(NrTS=NrTS, EndCriteria=1e-6)
FDTD.SetSinusExcite(f0)
FDTD.SetBoundaryCond(['PEC', 'PEC', 'PMC', 'PMC', 'PML_8', 'PML_8'])

CSX = ContinuousStructure()
FDTD.SetCSX(CSX)
mesh = CSX.GetGrid()
mesh.SetDeltaUnit(unit)
mesh.AddLine('x', [0, width/2, width])
mesh.AddLine('y', [0, width/2, width])
mesh.AddLine('z', np.arange(0, length + dz, dz))

CSX.AddExcitation('plane', exc_type=0, exc_val=[1, 0, 0]).AddBox([0, 0, z_src], [width, width, z_src])
CSX.AddProbe('E_probe', p_type=2).AddPoint([width/2, width/2, z_probe])

FDTD.Run(Sim_Path, cleanup=True)

### the timestep, from the excitation signal file: one sample per timestep
dT = np.loadtxt(os.path.join(Sim_Path, 'et'), comments='%')[1, 0]
TS = int(round(T0/dT))
print(f'  {TS} timesteps per period')
assert TS == TS0, f'FAIL: {TS} timesteps per period, expected {TS0} (the mesh or the timestep changed?)'

# field probe file columns: t/s, Ex, Ey, Ez
data = np.loadtxt(os.path.join(Sim_Path, 'E_probe'), comments='%')
t, val = data[:, 0], data[:, 1]
t_end = t[-1]
print(f'  simulation ended at {t_end*1e9:.2f} ns = {t_end/T0:.0f} periods')
assert t_end < 50*T0, \
    f'FAIL: simulation ran for {t_end/T0:.0f} periods, expected the steady-state detection to end it within 50'


def phasor(t_start, t_stop):
    """ least-squares fit of a cos(w t) + b sin(w t) in [t_start, t_stop), returns a - jb """
    m = (t >= t_start) & (t < t_stop)
    w = 2*np.pi*f0
    A = np.column_stack((np.cos(w*t[m]), np.sin(w*t[m])))
    (a, b), *_ = np.linalg.lstsq(A, val[m], rcond=None)
    return a - 1j*b


last = phasor(t_end - T0,   t_end)
prev = phasor(t_end - 2*T0, t_end - T0)
change = abs(last - prev) / abs(last)
print(f'  change between the last two periods: {change*100:.4f} %')
assert change < 1e-3, \
    f'FAIL: the probe changes {change*100:.3f} % from one period to the next, expected < 0.1 %'

print('PASS')
