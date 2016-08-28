# -*- coding: utf-8 -*-
"""
Created on Fri Dec 18 20:56:53 2015

@author: thorsten
"""

import os, tempfile
from pylab import *

from CSXCAD import CSXCAD

from openEMS.openEMS import openEMS
from openEMS.physical_constants import *

Sim_Path = os.path.join(tempfile.gettempdir(), 'Simp_Patch')

post_proc_only = True

# patch width in x-direction
patch_width  = 32 # resonant length
# patch length in y-direction
patch_length = 40

#substrate setup
substrate_epsR   = 3.38
substrate_kappa  = 1e-3 * 2*pi*2.45e9 * EPS0*substrate_epsR
substrate_width  = 60
substrate_length = 60
substrate_thickness = 1.524
substrate_cells = 4

#setup feeding
feed_pos = -6 #feeding position in x-direction
feed_R = 50     #feed resistance

# size of the simulation box
SimBox = np.array([200, 200, 150])

## setup FDTD parameter & excitation function
f0 = 2e9 # center frequency
fc = 1e9 # 20 dB corner frequency
FDTD = openEMS(NrTS=30000, EndCriteria=1e-4)
FDTD.SetGaussExcite( f0, fc )
FDTD.SetBoundaryCond( ['MUR', 'MUR', 'MUR', 'MUR', 'MUR', 'MUR'] )

mesh_res = C0/(f0+fc)/1e-3/20

CSX = CSXCAD.ContinuousStructure()
FDTD.SetCSX(CSX)


#initialize the mesh with the "air-box" dimensions
mesh = {}
mesh['x'] = [-SimBox[0]/2, feed_pos, SimBox[0]/2]
mesh['y'] = [-SimBox[1]/2, SimBox[1]/2]
mesh['z'] = [-SimBox[2]/3, SimBox[2]*2/3]

## create patch
patch = CSX.AddMetal( 'patch' ) # create a perfect electric conductor (PEC)
start = [-patch_width/2, -patch_length/2, substrate_thickness]
stop  = [ patch_width/2 , patch_length/2, substrate_thickness]
pb=CSX.AddBox(patch, priority=10, start=start, stop=stop) # add a box-primitive to the metal property 'patch'

edge_mesh = np.array([-1/3.0, 2/3.0]) * mesh_res/2
mesh['x'] = r_[mesh['x'], start[0]-edge_mesh, stop[0]+edge_mesh]
mesh['y'] = r_[mesh['y'], start[1]-edge_mesh, stop[1]+edge_mesh]

## create substrate
substrate = CSX.AddMaterial( 'substräte', Epsilon=substrate_epsR, Kappa=substrate_kappa)
start = [-substrate_width/2, -substrate_length/2, 0]
stop  = [ substrate_width/2,  substrate_length/2, substrate_thickness]
sb=CSX.AddBox( substrate, priority=0, start=start, stop=stop )

# add extra cells to discretize the substrate thickness
mesh['z'] = r_[mesh['z'], linspace(0,substrate_thickness,substrate_cells+1)]

## create ground (same size as substrate)
gnd = CSX.AddMetal( 'gnd' ) # create a perfect electric conductor (PEC)
start[2]=0
stop[2] =0
gb=CSX.AddBox(gnd, start, stop, priority=10)

mesh['x'] = r_[mesh['x'], start[0], stop[0]]
mesh['y'] = r_[mesh['y'], start[1], stop[1]]

## apply the excitation & resist as a current source
start = [feed_pos, 0, 0]
stop  = [feed_pos, 0, substrate_thickness]
port = FDTD.AddLumpedPort(1 ,feed_R, start, stop, 'z', 1.0, priority=5)

mesh['x'] = r_[mesh['x'], start[0]]
mesh['y'] = r_[mesh['y'], start[1]]

CSX.DefineGrid(mesh, unit=1e-3, smooth_mesh_res=mesh_res)

nf2ff = FDTD.CreateNF2FFBox()

if 0:  # debugging only
    CSX_file = os.path.join(Sim_Path, 'simp_patch.xml')
    CSX.Write2XML(CSX_file)
    os.system(r'AppCSXCAD "{}"'.format(CSX_file))

if not post_proc_only:
    FDTD.Run(Sim_Path, verbose=3, cleanup=True)

f = np.linspace(max(1e9,f0-fc),f0+fc,401)
port.CalcPort(Sim_Path, f)
s11 = port.uf_ref/port.uf_inc
s11_dB = 20.0*np.log10(np.abs(s11))
figure()
plot(f/1e9, s11_dB)
grid()
ylabel('s11 (dB)')
xlabel('frequency (GHz)')

idx = np.where((s11_dB<-10) & (s11_dB==np.min(s11_dB)))[0]
if not len(idx)==1:
    print('No resonance frequency found for far-field calulation')
else:
    f_res = f[idx[0]]
    theta = np.arange(-180.0, 180.0, 2.0)
    phi   = [0., 90.]
    nf2ff.CalcNF2FF(Sim_Path, f_res, theta, phi, center=[0,0,1e-3], read_cached=False )

    figure()
    E_norm = 20.0*np.log10(nf2ff.E_norm[0]/np.max(nf2ff.E_norm[0])) + nf2ff.Dmax[0]
    plot(theta, np.squeeze(E_norm[:,0]), label='xz-plane')
    plot(theta, np.squeeze(E_norm[:,1]), label='yz-plane')
    grid()
    ylabel('directivity (dBi)')
    xlabel('theta (deg)')
    title('Frequency: {} GHz'.format(f_res/1e9))
    legend()

Zin = port.uf_tot/port.if_tot
figure()
plot(f/1e9, np.real(Zin))
plot(f/1e9, np.imag(Zin))
grid()
ylabel('Zin (Ohm)')
xlabel('frequency (GHz)')

show()
