"""
 Microstrip Line with local absorber test

 (c) 20@3-2025 Gadi Lahav <gadi@rfwithcare.com>

"""


### Import Libraries
import os, tempfile
from pylab import *
import scipy.io

from CSXCAD  import ContinuousStructure
from openEMS import openEMS
from openEMS.physical_constants import *

from CSXCAD.CSProperties import ABCtype

### General parameter setup
Sim_Path = os.path.join(tempfile.gettempdir(), 'Test_MSL_W_SA')

print(Sim_Path)

post_proc_only = False

# patch width (resonant length) in x-direction
microstrip_W = 1.875 #

#substrate setup
substrate_epsR   = 4.5
substrate_width  = 15
substrate_length = 50
substrate_thickness = 1
substrate_cells = 8

mesh_res = 0.1
# mesh_res = 5

# Port setup
port_w_fact = 6
port_h_fact = 5.5
port_thick_mm = 0.5
port_shift_mm = 15


cu_thick = 0.1

Airbox_Add = 12.5;


# setup FDTD parameter & excitation function
f0 = 2e9 # center frequency
fc = 1e9 # 20 dB corner frequency


# size of the simulation box
SimBox = np.array([
            -substrate_width*0.5 - Airbox_Add,
            substrate_width*0.5 + Airbox_Add,
            -cu_thick - Airbox_Add, 
            substrate_thickness*(1 + port_h_fact) + Airbox_Add,
            -Airbox_Add,
            substrate_length + Airbox_Add])


### FDTD setup
## * Limit the simulation to 30k timesteps
## * Define a reduced end criteria of -40dB
FDTD = openEMS(NrTS=30000, EndCriteria=1e-4)
FDTD.SetGaussExcite( f0, fc )
# FDTD.SetBoundaryCond( ['MUR', 'MUR', 'MUR', 'MUR', 'MUR', 'MUR'] )
FDTD.SetBoundaryCond( ['PML_8', 'PML_8','PML_8','PML_8','PML_8','PML_8']);

CSX = ContinuousStructure()
FDTD.SetCSX(CSX)
mesh = CSX.GetGrid()
mesh.SetDeltaUnit(1e-3)
mesh_res = ((C0/(f0+fc))/1e-3)/75

### Generate properties, primitives and mesh-grid
#initialize the mesh with the "air-box" dimensions
mesh.AddLine('x', SimBox[0:2])
mesh.AddLine('y', SimBox[2:4])
mesh.AddLine('z', SimBox[4:6])

# create line and ground
line = CSX.AddMetal('cu_top')
start = [-microstrip_W/2, substrate_thickness, 0.0]
stop  = [ microstrip_W/2, substrate_thickness + cu_thick, substrate_length]
line.AddBox(priority = 10, start = start, stop = stop) # add a box-primitive to the metal property 'patch'
FDTD.AddEdges2Grid(dirs='xyz', properties=line)
mesh.AddLine('x',[start[0],stop[0]])
mesh.AddLine('y',[start[1],stop[1]])
mesh.AddLine('z',[start[2],stop[2]])

# Port PEC block
port2_block = CSX.AddMetal('PEC')
start = [-microstrip_W*(1.0 + port_w_fact)*0.5, -cu_thick, substrate_length - port_shift_mm]
stop  = [ microstrip_W*(1.0 + port_w_fact)*0.5, substrate_thickness*(1.0 + port_h_fact), substrate_length + port_thick_mm - port_shift_mm]
port2_block.AddBox(priority=12, start=start, stop=stop) # add a box-primitive to the metal property 'patch'
FDTD.AddEdges2Grid(dirs='xyz', properties=port2_block)
mesh.AddLine('x',[start[0],stop[0]])
mesh.AddLine('y',[start[1],stop[1]])
mesh.AddLine('z',[start[2],stop[2]])

# Extra pad for port not to take a lot of space
mesh.AddLine('z',substrate_length - port_shift_mm - np.array([0.5,1])*mesh_res)

# create substrate
sub = CSX.AddMaterial('FR4', epsilon = substrate_epsR)
start = [-substrate_width/2, 0.0, 0.0]
stop  = [ substrate_width/2, substrate_thickness,  substrate_length]
sub.AddBox(priority = 2, start = start, stop = stop)
mesh.AddLine('x',[start[0],stop[0]])
mesh.AddLine('y',[start[1],stop[1]])
mesh.AddLine('z',[start[2],stop[2]])


# add extra cells to discretize the substrate thickness
mesh.AddLine('y', linspace(0, substrate_thickness, substrate_cells+1))

# create ground (same size as substrate)
gnd = CSX.AddMaterial('cu_bot',kappa=56000000)
start = [-substrate_width/2, -cu_thick, 0.0]
stop  = [ substrate_width/2, 0.0, substrate_length]
gnd.AddBox(priority=10, start=start, stop=stop) # add a box-primitive to the metal property 'patch'
FDTD.AddEdges2Grid(dirs='xyz', properties=gnd)

mesh.SmoothMeshLines('all', mesh_res, 1.4)

# Find start\stop mesh lines
Zz = mesh.GetLines('z')
idxTerm = (np.where(Zz == (substrate_length - port_shift_mm))[0]).item(0)

# apply the excitation & resist as a current source
# Define port mode 
start = [-microstrip_W*0.5, 0.0, 0.0]
stop  = [ microstrip_W*0.5, substrate_thickness, 0.0]
port1 = FDTD.AddLumpedPort(1, 50.0, start, stop, 'y', 1.0, priority=15, edges2grid='xy')


# Place absorber for port #2

v_phase = np.sqrt(1/(0.5*(1 + substrate_epsR)))*C0

# Use this option to place directly on PEC. Absorption is ~ 3dB worse, this way, but less leakage.
start = [-microstrip_W*(1.0 + port_w_fact)*0.5, -cu_thick, Zz.item(idxTerm)]
stop  = [ microstrip_W*(1.0 + port_w_fact)*0.5, substrate_thickness*(1.0 + port_h_fact), Zz.item(idxTerm)]
# This places the absorber one mesh cell farther. 
# start = [-microstrip_W*(1.0 + port_w_fact)*0.5, -cu_thick, Zz.item(idxTerm - 1)]
# stop  = [ microstrip_W*(1.0 + port_w_fact)*0.5, substrate_thickness*(1.0 + port_h_fact), Zz.item(idxTerm - 1)]

# In this application, MUR_1ST_SA may cause instability
abs2 = CSX.AddAbsorbingBC('abs2',NormalSignPositive = False, AbsorbingBoundaryType = ABCtype.MUR_1ST)
abs2.AddBox(start, stop, priority=20)

abs2.SetPhaseVelocity(v_phase);

# # Define dump box...W
# Et = CSX.AddDump('Et', file_type=0, dump_type=0, dump_mode=1)
# start = [float(SimBox[0]), float(SimBox[2]), float(SimBox[4])];
# stop  = [float(SimBox[1]), float(SimBox[3]), float(SimBox[5])];
# Et.AddBox(start, stop);


### Run the simulation
if 1:  # debugging only
    CSX_file = os.path.join(Sim_Path, 'msl_w_lumped_port_and_absorber.xml')
    if not os.path.exists(Sim_Path):
        os.mkdir(Sim_Path)
    CSX.Write2XML(CSX_file)
    from CSXCAD import AppCSXCAD_BIN
    os.system(AppCSXCAD_BIN + ' "{}"'.format(CSX_file))
    

if not post_proc_only:
    FDTD.Run(Sim_Path, verbose=0, cleanup=False)

### Post-processing and plrorotting
f = np.linspace(max(1e9,f0-fc),f0+fc,401)
port1.CalcPort(Sim_Path, f, ref_impedance = 50.0)
s11 = port1.uf_ref/port1.uf_inc

s11_dB = np.transpose(20.0*np.log10(np.abs(s11)))

figure()
plot(f/1e9, s11_dB, 'k-', linewidth=2, label='$S_{11}$')
grid()
legend()
ylabel('S-Parameter (dB)')
xlabel('Frequency (GHz)')

show()

