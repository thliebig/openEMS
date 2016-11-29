# -*- coding: utf-8 -*-
"""
Created on Tue Dec 29 15:40:23 2015

 Tutorials / Rect_Waveguide

 Describtion at:
 http://openems.de/index.php/Tutorial:_Rectangular_Waveguide

 Tested with
  - python 3.4
  - openEMS v0.0.33+

 (C) 2015-2016 Thorsten Liebig <thorsten.liebig@gmx.de>

"""

import os, tempfile
from pylab import *

from CSXCAD import CSXCAD
from openEMS.openEMS import openEMS
from openEMS.physical_constants import *

Sim_Path = os.path.join(tempfile.gettempdir(), 'Rect_WG')

## setup the simulation ###################################################
post_proc_only = False
unit = 1e-6; #drawing unit in um

# waveguide dimensions
# WR42
a = 10700;   #waveguide width
b = 4300;    #waveguide heigth
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

## setup FDTD parameter & excitation function #############################
FDTD = openEMS(NrTS=1e4);
FDTD.SetGaussExcite(0.5*(f_start+f_stop),0.5*(f_stop-f_start));

# boundary conditions
FDTD.SetBoundaryCond([0, 0, 0, 0, 3, 3]);

## setup CSXCAD geometry & mesh
CSX = CSXCAD.ContinuousStructure()
FDTD.SetCSX(CSX)
mesh = CSX.GetGrid()
mesh.SetDeltaUnit(unit)

mesh.AddLine('x', [0, a])
mesh.AddLine('y', [0, b])
mesh.AddLine('z', [0, length])

## apply the waveguide port ###################################################
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

## define dump box... #####################################################
Et = CSX.AddDump('Et', file_type=0, sub_sampling=[2,2,2])
start = [0, 0, 0];
stop  = [a, b, length];
CSX.AddBox(Et, start, stop);

## Write openEMS compatoble xml-file ######################################
if 0:  # debugging only
    CSX_file = os.path.join(Sim_Path, 'rect_wg.xml')
    CSX.Write2XML(CSX_file)
    os.system(r'AppCSXCAD "{}"'.format(CSX_file))

if not post_proc_only:
    FDTD.Run(Sim_Path, verbose=3, cleanup=True)

### postproc ###############################################################
freq = linspace(f_start,f_stop,201)
for port in ports:
    port.CalcPort(Sim_Path, freq)

s11 = ports[0].uf_ref / ports[0].uf_inc
s21 = ports[1].uf_ref / ports[0].uf_inc
ZL  = ports[0].uf_tot / ports[0].if_tot
ZL_a = ports[0].ZL # analytic waveguide impedance

## plot s-parameter #######################################################
figure()
plot(freq*1e-6,20*log10(abs(s11)),'k-',linewidth=2, label='$S_{11}$')
grid()
plot(freq*1e-6,20*log10(abs(s21)),'r--',linewidth=2, label='$S_{21}$')
legend();
ylabel('S-Parameter (dB)')
xlabel(r'frequency (MHz) $\rightarrow$')

## compare analytic and numerical wave-impedance ##########################
figure()
plot(freq*1e-6,real(ZL), linewidth=2, label='$\Re\{Z_L\}$')
grid()
plot(freq*1e-6,imag(ZL),'r--', linewidth=2, label='$\Im\{Z_L\}$')
plot(freq*1e-6,ZL_a,'g-.',linewidth=2, label='$Z_{L, analytic}$')
ylabel('ZL $(\Omega)$')
xlabel(r'frequency (MHz) $\rightarrow$')
legend()

show()
