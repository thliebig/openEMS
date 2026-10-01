# -*- coding: utf-8 -*-
"""
 Rectangular Waveguide with local absorbers Test

 A WR42 section closed by PEC at both ends, with a local Mur absorbing boundary
 (MUR_1ST_SA) placed in front of each end. Verifies that the local absorber
 terminates the guide -- the PEC box would be a resonator without it -- and that
 the wave impedance of AddRectWaveGuidePort matches the analytic TE10 value.

 Pass criteria:
   max(dB(S11)) < -40 dB    (the absorbers terminate the guide; measured -48 dB,
                             a reflecting end would give ~0 dB)
   |dB(S21)| < 0.05 dB      (lossless section between the two ports)
   real(ZL) within 3 % of the analytic TE10 wave impedance, |imag(ZL)| below
   3 % of it (measured 0.8 % each at this lambda/50 mesh)

 Tested with
  - python 3.14
  - openEMS v0.0.36+

 (c) 2023-2025 Gadi Lahav <gadi@rfwithcare.com>

"""

### Import Libraries
import os, tempfile
import numpy as np

from CSXCAD  import ContinuousStructure
from openEMS import openEMS
from openEMS.physical_constants import *

from CSXCAD.CSProperties import ABCtype

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
# mesh_res = lambda0/30
mesh_res = lambda0/50

### Setup FDTD parameter & excitation function
FDTD = openEMS(NrTS=1e4);
FDTD.SetGaussExcite(0.5*(f_start+f_stop),0.5*(f_stop-f_start));

# boundary conditions
FDTD.SetBoundaryCond([0, 0, 0, 0, 0, 0]);

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
start=[0, 0, 3*mesh_res];
stop =[a, b, 4*mesh_res];
mesh.AddLine('z', [start[2], stop[2]])
ports.append(FDTD.AddRectWaveGuidePort( 0, start, stop, 'z', a*unit, b*unit, TE_mode, 1))

start=[0, 0, length-3*mesh_res];
stop =[a, b, length-4*mesh_res];
mesh.AddLine('z', [start[2], stop[2]])
ports.append(FDTD.AddRectWaveGuidePort( 1, start, stop, 'z', a*unit, b*unit, TE_mode))

# Add PEC Boxes
pecBlocks = CSX.AddMetal('PEC')

start = [0, 0, 0.0]
stop  = [ a, b, 1*mesh_res]
pecBlocks.AddBox(priority=5, start=start, stop=stop) # add a box-primitive to the metal property 'patch'
start = [ 0, 0, length-1*mesh_res]
stop  = [a, b, length]
pecBlocks.AddBox(priority=5, start=start, stop=stop) # add a box-primitive to the metal property 'patch'

start = [0, 0., 2*mesh_res]
stop  = [a, b, 2*mesh_res]

mesh.AddLine('z', [stop[2]])
abs1 = CSX.AddAbsorbingBC('abs1',NormalSignPositive = True, AbsorbingBoundaryType = ABCtype.MUR_1ST_SA, PhaseVelocity = 3.6196e+08)
abs1.AddBox(start, stop, priority=6)



start=[0, 0, length-2*mesh_res];
stop =[a, b, length-2*mesh_res];

mesh.AddLine('z', [stop[2]])
abs2 = CSX.AddAbsorbingBC('abs2',NormalSignPositive = False, AbsorbingBoundaryType = ABCtype.MUR_1ST_SA, PhaseVelocity = 3.6196e+08)
abs2.AddBox(start, stop, priority=6)

mesh.SmoothMeshLines('all', mesh_res, ratio=1.4)

### Define dump box...
# Et = CSX.AddDump('Et', file_type=0, sub_sampling=[2,2,2])
# start = [0, 0, 0];
# stop  = [a, b, length];
# Et.AddBox(start, stop);

### Run the simulation
if 0:  # set to 1 to inspect the geometry in AppCSXCAD
    CSX_file = os.path.join(Sim_Path, 'rect_wg.xml')
    if not os.path.exists(Sim_Path):
        os.mkdir(Sim_Path)
    CSX.Write2XML(CSX_file)
    from CSXCAD import AppCSXCAD_BIN
    os.system(AppCSXCAD_BIN + ' "{}"'.format(CSX_file))

if not post_proc_only:
    FDTD.Run(Sim_Path, cleanup=True, exact_endcriteria=True)

### Post-processing
freq = np.linspace(f_start, f_stop, 201)
for port in ports:
    port.CalcPort(Sim_Path, freq)

s11 = ports[0].uf_ref / ports[0].uf_inc
s21 = ports[1].uf_ref / ports[0].uf_inc
ZL   = ports[0].uf_tot / ports[0].if_tot
ZL_a = ports[0].ZL  # analytic TE10 wave impedance

s11_dB = 20*np.log10(np.abs(s11))
s21_dB = 20*np.log10(np.abs(s21))

### Pass / fail checks
print('max(dB(S11))  = {:.1f} dB'.format(np.max(s11_dB)))
print('dB(S21)       = {:.3f} .. {:.3f} dB'.format(np.min(s21_dB), np.max(s21_dB)))

# Both ends are PEC with a local absorber in front of them, so the only way out
# of the box is through the absorbers: a broken one turns this into a resonator.
assert np.max(s11_dB) < -40, \
    'FAIL: max(dB(S11)) = {:.1f} dB, expected < -40 dB'.format(np.max(s11_dB))

# lossless section between the ports -- and a guard against a sign or
# normalisation error, which would show up as gain
assert np.min(s21_dB) > -0.05, \
    'FAIL: min(dB(S21)) = {:.3f} dB, expected > -0.05 dB'.format(np.min(s21_dB))
assert np.max(s21_dB) < 0.05, \
    'FAIL: max(dB(S21)) = {:.3f} dB, expected < +0.05 dB (sign error?)'.format(np.max(s21_dB))

# the wave impedance the port reports against the closed form of the TE10 mode;
# the tolerance is what the lambda/50 mesh and the de-embedding leave over
err_re = np.max(np.abs(np.real(ZL) - ZL_a) / ZL_a)
err_im = np.max(np.abs(np.imag(ZL)) / ZL_a)
print('ZL error      = {:.1f} % real, {:.1f} % imaginary'.format(err_re*100, err_im*100))
assert err_re < 0.03, \
    'FAIL: real(ZL) deviates by {:.1f} %, expected < 3 %'.format(err_re*100)
assert err_im < 0.03, \
    'FAIL: imag(ZL) reaches {:.1f} % of ZL, expected < 3 %'.format(err_im*100)

print('PASS')

if 0:  # set to 1 for debugging plots
    import matplotlib.pyplot as plt

    fig, axis = plt.subplots(num='S-Parameters', tight_layout=True)
    axis.plot(freq/1e9, s11_dB, 'k-',  linewidth=2, label='$S_{11}$')
    axis.plot(freq/1e9, s21_dB, 'r--', linewidth=2, label='$S_{21}$')
    axis.grid()
    axis.set_xmargin(0)
    axis.set_xlabel('Frequency (GHz)')
    axis.set_ylabel('S-Parameter (dB)')
    axis.legend()

    fig, axis = plt.subplots(num='Wave Impedance', tight_layout=True)
    axis.plot(freq/1e9, np.real(ZL), linewidth=2, label=r'$\Re\{Z_L\}$')
    axis.plot(freq/1e9, np.imag(ZL), 'r--', linewidth=2, label=r'$\Im\{Z_L\}$')
    axis.plot(freq/1e9, ZL_a, 'g-.', linewidth=2, label='$Z_{L,analytic}$')
    axis.grid()
    axis.set_xmargin(0)
    axis.set_xlabel('Frequency (GHz)')
    axis.set_ylabel(r'Wave impedance $(\Omega)$')
    axis.legend()

    plt.show()
