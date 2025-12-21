# -*- coding: utf-8 -*-
"""
 Rectangular Waveguide Tutorial

 Tested with
  - python 3.10
  - openEMS v0.0.35+

 (c) 2015-2023 Thorsten Liebig <thorsten.liebig@gmx.de>
 15-Dec-2025: modified to use matplotlib.pyplot instead of pylab

"""

### Import Libraries
import os, tempfile
import numpy as np
import matplotlib.pyplot as plt  # pip install matplotlib

from CSXCAD  import ContinuousStructure
from openEMS import openEMS
from openEMS.physical_constants import *

### Setup the simulation
Sim_Path = os.path.join(tempfile.gettempdir(), 'Rect_WG')

post_proc_only = False
unit = 1e-6; #drawing unit in um

# waveguide dimensions
# WR42
a = 10700;   #waveguide width
b = 4300;    #waveguide height
length = 50000;

# frequency range of interest
f_start = 20e9;
f_0     = 24e9;
f_stop  = 26e9;
lambda0 = C0/f_0/unit;

#waveguide TE-mode definition
TE_mode = 'TE10';

#targeted mesh resolution
mesh_res = lambda0/30

### Setup FDTD parameter & excitation function
FDTD = openEMS(NrTS=1e4);
FDTD.SetGaussExcite(0.5*(f_start+f_stop),0.5*(f_stop-f_start));

# boundary conditions
FDTD.SetBoundaryCond([0, 0, 0, 0, 3, 3]);

### Setup geometry & mesh
CSX = ContinuousStructure()
FDTD.SetCSX(CSX)
mesh = CSX.GetGrid()
mesh.SetDeltaUnit(unit)

mesh.AddLine('x', [0, a])
mesh.AddLine('y', [0, b])
mesh.AddLine('z', [0, length])

## Apply the waveguide port
ports = []
start=[0, 0, 10*mesh_res];
stop =[a, b, 15*mesh_res];
mesh.AddLine('z', [start[2], stop[2]])
ports.append(FDTD.AddRectWaveGuidePort( 0, start, stop, 'z', a*unit, b*unit, TE_mode, 1))

start=[0, 0, length-10*mesh_res];
stop =[a, b, length-15*mesh_res];
mesh.AddLine('z', [start[2], stop[2]])
ports.append(FDTD.AddRectWaveGuidePort( 1, start, stop, 'z', a*unit, b*unit, TE_mode))

mesh.SmoothMeshLines('all', mesh_res, ratio=1.4)

### Define dump box...
Et = CSX.AddDump('Et', file_type=0, sub_sampling=[2,2,2])
start = [0, 0, 0];
stop  = [a, b, length];
Et.AddBox(start, stop);

### Run the simulation
if 0:  # debugging only
    CSX_file = os.path.join(Sim_Path, 'rect_wg.xml')
    if not os.path.exists(Sim_Path):
        os.mkdir(Sim_Path)
    CSX.Write2XML(CSX_file)
    from CSXCAD import AppCSXCAD_BIN
    os.system(AppCSXCAD_BIN + ' "{}"'.format(CSX_file))

if not post_proc_only:
    FDTD.Run(Sim_Path, cleanup=True)

### Postprocessing & plotting
freq = np.linspace(f_start,f_stop,201)
for port in ports:
    port.CalcPort(Sim_Path, freq)

s11 = ports[0].uf_ref / ports[0].uf_inc
s21 = ports[1].uf_ref / ports[0].uf_inc
ZL  = ports[0].uf_tot / ports[0].if_tot
ZL_a = ports[0].ZL # analytic waveguide impedance

## Plot s-parameter

fig, axis = plt.subplots(num="S11", tight_layout=True)
axis.plot(freq/1e9, 20*np.log10(abs(s11)), 'k-',  linewidth=2, label='dB(S11)')
axis.plot(freq/1e9, 20*np.log10(abs(s21)), 'r--',  linewidth=2, label='dB(S21)')
axis.grid()
axis.set_xmargin(0)
axis.set_xlabel('Frequency (GHz)')
axis.set_ylabel('S-Parameter (dB)')
axis.legend()


## Compare analytic and numerical wave-impedance
fig, axis = plt.subplots(num="ZL", tight_layout=True)
axis.plot(freq/1e9, np.real(ZL), 'k-',  linewidth=2, label='Re(ZL)')
axis.plot(freq/1e9, np.imag(ZL), 'r--', linewidth=2, label='Im(ZL)')

axis.plot(freq/1e9, ZL_a, 'g-.',  linewidth=2, label='Re(ZL, analytic)')
axis.grid()
axis.set_xmargin(0)
axis.set_xlabel('Frequency (GHz)')
axis.set_ylabel('ZL (Ohm)')
axis.legend()

# show all plots
plt.show()
