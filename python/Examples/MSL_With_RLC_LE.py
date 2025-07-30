######################
######################
## Import Libraries ##
######################
######################

import os, tempfile
from pylab import *

from CSXCAD  import ContinuousStructure
from openEMS import openEMS
from openEMS.physical_constants import *
from CSXCAD.CSProperties import LEtype
from cmath import pi

import skrf as rf

########################
########################
## Project Parameters ##
########################
########################

# Only post processing
post_proc_only = False

# PCB settings
substrate_epsR   = 4.5
substrate_width  = 8
substrate_length = 10
substrate_thickness = 1
substrate_cells = 5
cu_thick = 0.1
microstrip_W = 1.875 

# Add airbox gap, in mm
Airbox_Add = 7.5;

# RLC lumped-element values
R = 10
L = 2e-9
C = 1e-12
Z0 = 50

networkType = LEtype.LE_SERIES # Can be LEtype.LE_SERIES or LEtype.LE_PARALLEL
 
# setup FDTD parameter & excitation function
f0 = 2e9 # center frequency
fc = 1e9 # 20 dB corner frequency

#################
#################
## Model Setup ##
#################
#################

### General parameter setup
Sim_Path = os.path.join(tempfile.gettempdir(), 'Simp_MSL_W_RLC')

# size of the simulation box
SimBox = np.array([
            -substrate_width*0.5 - Airbox_Add,
            substrate_width*0.5 + Airbox_Add,
            -Airbox_Add,
            substrate_length + Airbox_Add, 
            -cu_thick - Airbox_Add, 
            substrate_thickness + cu_thick + Airbox_Add])



### FDTD setup
## * Limit the simulation to 30k timesteps
## * Define a reduced end criteria of -40dB
FDTD = openEMS(NrTS=160000, EndCriteria=1e-4)
FDTD.SetGaussExcite( f0, fc )
FDTD.SetBoundaryCond( ['MUR', 'MUR', 'MUR', 'MUR', 'MUR', 'MUR'] )


CSX = ContinuousStructure()
FDTD.SetCSX(CSX)
mesh = CSX.GetGrid()
mesh.SetDeltaUnit(1e-3)
mesh_res = C0/(f0+fc)/1e-3/150

### Generate properties, primitives and mesh-grid
#initialize the mesh with the "air-box" dimensions
mesh.AddLine('x', SimBox[0:2])
mesh.AddLine('y', SimBox[2:4])
mesh.AddLine('z', SimBox[4:6])

# create line and ground
line = CSX.AddMaterial('cu_top',kappa=56000000)
start = [-microstrip_W/2, 0.0, substrate_thickness]
stop  = [ microstrip_W/2 , substrate_length, substrate_thickness + cu_thick]
line.AddBox(priority=10, start=start, stop=stop) # add a box-primitive to the metal property 'patch'
FDTD.AddEdges2Grid(dirs='xyz', properties=line)

# create substrate
sub = CSX.AddMaterial( 'FR4', epsilon=substrate_epsR)
start = [-substrate_width/2, 0.0, 0.0]
stop  = [ substrate_width/2,  substrate_length, substrate_thickness]
sub.AddBox( priority=2, start=start, stop=stop )

# add extra cells to discretize the substrate thickness
mesh.AddLine('z', linspace(0,substrate_thickness,substrate_cells+1))

# create ground (same size as substrate)
gnd = CSX.AddMaterial('cu_bot',kappa=56000000)
start = [-substrate_width/2, 0.0, -cu_thick]
stop  = [ substrate_width/2, substrate_length, 0.0]
gnd.AddBox(priority=10, start=start, stop=stop) # add a box-primitive to the metal property 'patch'
FDTD.AddEdges2Grid(dirs='xyz', properties=gnd)

# apply the excitation & resist as a current source
start = [-microstrip_W*0.5, 0, 0]
stop  = [ microstrip_W*0.5, 0, substrate_thickness]
port = FDTD.AddLumpedPort(1, Z0, start, stop, 'z', 1.0, priority=15, edges2grid='xy')

start = [-microstrip_W*0.5, substrate_length, 0]
stop  = [ microstrip_W*0.5, substrate_length, substrate_thickness]
LE = CSX.AddLumpedElement('RLC_SER', ny='z', caps=False, R=R, L=L, C=C,LEtype=networkType)
LE.AddBox(start, stop, priority=25)


mesh.SmoothMeshLines('all', mesh_res, 1.4)


### Run the simulation
if 0:  # debugging only
    CSX_file = os.path.join(Sim_Path, 'simp_msl_w_le.xml')
    if not os.path.exists(Sim_Path):
        os.mkdir(Sim_Path)
    CSX.Write2XML(CSX_file)
    from CSXCAD import AppCSXCAD_BIN
    os.system(AppCSXCAD_BIN + ' "{}"'.format(CSX_file))

if not post_proc_only:
    FDTD.Run(Sim_Path, verbose=3, cleanup=True)

##################################
##################################
## Post-processing and plotting ##
##################################
##################################

f = np.linspace(max(1e9,f0-fc),f0+fc,401)
port.CalcPort(Sim_Path, f)
s11 = port.uf_ref/port.uf_inc
s11_dB = 20.0*np.log10(np.abs(s11))

# Analytical formulation of RLC network
if networkType == LEtype.LE_SERIES:
    if not (C == 0):
        Zref  = R + 1j*2*pi*f*L + 1.0/(1j*2.0*pi*f*C)
    else:
        Zref  = R + 1j*2*pi*f*L
elif networkType == LEtype.LE_PARALLEL:
    if not (L == 0):
        Yref  = 1/R + 1/(1j*2*pi*f*L) + 1j*2.0*pi*f*C
    else:
        Yref  = 1/R + 1j*2.0*pi*f*C
    Zref = 1/Yref
else:
    sys.exit("Lumped-Element network type undefined")

Gref = (Zref - Z0)/(Zref + Z0)
Sref = Gref.reshape(-1, 1, 1)

# De-embed simulation, using naive matched line de-embedding
Dk_eff = 0.5*(substrate_epsR + 1.0)
s11 = s11/exp(-2*1j*substrate_length*1e-3*2*pi*f/(0.87*C0/sqrt(Dk_eff)))

# Setup smith plot
if networkType == LEtype.LE_SERIES:
    networkName = "Series RLC - Analytical"
else:
    networkName = "Parallel RLC - Analytical"
    
ntw = rf.Network(frequency=f/1e9, s=Sref, z0=50, name=networkName)
plt.figure(figsize=(8, 6))
ntw.plot_s_smith()

# Add reference plot
plt.plot(s11.real, s11.imag, '--', label='openEMS - De-embedded, S11', color='red')

if networkType == LEtype.LE_SERIES:
    plt.title('Series RLC Reflection Coefficient')
else:
    plt.title('Parallel RLC Reflection Coefficient')
    
plt.legend()
plt.grid(True)
plt.show()
