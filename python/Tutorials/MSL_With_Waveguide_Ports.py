### Import Libraries
import os, tempfile
from pylab import *
import scipy.io

from CSXCAD  import ContinuousStructure
from openEMS import openEMS
from openEMS.physical_constants import *
from CSXCAD.CSProperties import BCtype

### General parameter setup
Sim_Path = os.path.join(tempfile.gettempdir(), 'Test_MSL_TEM_Solver')

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
portThick_mm = 0.5

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
FDTD.SetBoundaryCond( ['PML_8', 'PML_8', 'PML_8', 'PML_8', 'PML_8', 'PML_8'] )

CSX = ContinuousStructure()
FDTD.SetCSX(CSX)
mesh = CSX.GetGrid()
mesh.SetDeltaUnit(1e-3)
mesh_res = ((C0/(f0+fc))/1e-3)/100

### Generate propertiesprimitives and mesh-grid
#initialize the mesh with the "air-box" dimensions
mesh.AddLine('x', SimBox[0:2])
mesh.AddLine('y', SimBox[2:4])
mesh.AddLine('z', SimBox[4:6])

# create line and ground
# line = CSX.AddMaterial('cu_top',kappa=56000000)
line = CSX.AddMetal('cu_top')
start = [-microstrip_W/2, substrate_thickness, 0.0]
stop  = [ microstrip_W/2, substrate_thickness + cu_thick, substrate_length]
line.AddBox(priority = 10, start = start, stop = stop) # add a box-primitive to the metal property 'patch'
FDTD.AddEdges2Grid(dirs='xyz', properties=line)
mesh.AddLine('x',[start[0],stop[0]])
mesh.AddLine('y',[start[1],stop[1]])
mesh.AddLine('z',[start[2],stop[2]])

# PEC block for Z = 0
port1_block = CSX.AddMetal('PEC')
start = [-microstrip_W*(1.0 + port_w_fact)*0.5, 0.0, -portThick_mm]
stop  = [ microstrip_W*(1.0 + port_w_fact)*0.5, substrate_thickness*(1.0 + port_h_fact), 0.0]
# port1_block.AddBox(priority=12, start=start, stop=stop) # add a box-primitive to the metal property 'patch'
# FDTD.AddEdges2Grid(dirs='xyz', properties=port1_block)
mesh.AddLine('x',[start[0],stop[0]])
mesh.AddLine('y',[start[1],stop[1]])
mesh.AddLine('z',[start[2],stop[2]])

# Extra pad for port not to take a lot of space
mesh.AddLine('z',np.array([0.25,0.5,0.8,1])*mesh_res)

port2_block = CSX.AddMetal('PEC')
start = [-microstrip_W*(1.0 + port_w_fact)*0.5, 0.0, substrate_length]
stop  = [ microstrip_W*(1.0 + port_w_fact)*0.5, substrate_thickness*(1.0 + port_h_fact), substrate_length + portThick_mm]
# port2_block.AddBox(priority=12, start=start, stop=stop) # add a box-primitive to the metal property 'patch'
# FDTD.AddEdges2Grid(dirs='xyz', properties=port2_block)
mesh.AddLine('x',[start[0],stop[0]])
mesh.AddLine('y',[start[1],stop[1]])
mesh.AddLine('z',[start[2],stop[2]])


# Extra pad for port not to take a lot of space
mesh.AddLine('z',substrate_length - np.array([0.25,0.5,0.8,1])*mesh_res)

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
idxPort1 = (np.where(Zz == 0.0)[0] + 2).item(0)
idxPort2 = (np.where(Zz == substrate_length)[0] - 2).item(0)
idxAbs1 = idxPort1 - 1
idxAbs2 = idxPort2 + 1

matFile = scipy.io.loadmat("MSL_Mode_1.mat")
Em1 = matFile["E"]
Hm1 = matFile["H"]

# Define port mode 
start = [-microstrip_W*(1.0 + port_w_fact)*0.5, 0.0, Zz.item(idxPort1 + 0)]
stop  = [ microstrip_W*(1.0 + port_w_fact)*0.5, substrate_thickness*(1.0 + port_h_fact), Zz.item(idxPort1 + 1)]
port1 = FDTD.AddWaveGuidePort( 1, start, stop, 'z', Em1, Hm1, 0.0, excite=1, excite_type = 0)

# Place absorber for port #1
stop[2] = Zz.item(idxAbs1)
start[2] = stop[2] # Place absorber at end of port
abs1 = CSX.AddAbsorbingBC('abs1',NormDir = 3, BCtype = BCtype.MUR_1ST_1PV_SA)
abs1.AddBox(start, stop, priority = 20)

matFile = scipy.io.loadmat("MSL_Mode_2.mat")
Em2 = matFile["E"]
Hm2 = matFile["H"]
kz = matFile["Kz"]


# Place absorber for port #2
start[2] = Zz.item(idxAbs2)
stop[2] = start[2]
abs2 = CSX.AddAbsorbingBC('abs2',NormDir = -3, BCtype = BCtype.MUR_1ST_1PV_SA)
abs2.AddBox(start, stop, priority=20)

start = [-microstrip_W*(1.0 + port_w_fact)*0.5, 0.0, Zz.item(idxPort2 - 0)]
stop  = [ microstrip_W*(1.0 + port_w_fact)*0.5, substrate_thickness*(1.0 + port_h_fact), Zz.item(idxPort2 - 1)]
port2 = FDTD.AddWaveGuidePort( 2, start, stop, 'z', Em2, Hm2, 0.0, excite=0, excite_type = 0)

v_phase = 2*np.pi*f0/kz
v_phase = v_phase[0]
abs1.SetPhaseVelocity(v_phase);
abs2.SetPhaseVelocity(v_phase);

# # Define dump box...
# Et = CSX.AddDump('Et', file_type=0, dump_type=0, dump_mode=1)
# start = [float(SimBox[0]), float(SimBox[2]), float(SimBox[4])];
# stop  = [float(SimBox[1]), float(SimBox[3]), float(SimBox[5])];
# Et.AddBox(start, stop);


### Run the simulation
if 1:  # debugging only
    CSX_file = os.path.join(Sim_Path, 'simp_msl_2_WG_ports.xml')
    if not os.path.exists(Sim_Path):
        os.mkdir(Sim_Path)
    CSX.Write2XML(CSX_file)
    from CSXCAD import AppCSXCAD_BIN
    os.system(AppCSXCAD_BIN + ' "{}"'.format(CSX_file))

if not post_proc_only:
    FDTD.Run(Sim_Path, verbose=0, cleanup=False)


Zl = matFile["Zl"]
Zw = matFile["Zw"]

### Post-processing and plotting
f = np.linspace(max(1e9,f0-fc),f0+fc,401)
port1.CalcPort(Sim_Path, f, ref_impedance = Zw, ZL = Zl)
port2.CalcPort(Sim_Path, f, ref_impedance = Zw, ZL = Zl)
s11 = port1.uf_ref/port1.uf_inc
s21 = port2.uf_ref/port1.uf_inc

s11_dB = np.transpose(20.0*np.log10(np.abs(s11)))
s21_dB = np.transpose(20.0*np.log10(np.abs(s21)))

figure()
plot(f/1e9, s11_dB, 'k-', linewidth=2, label='$S_{11}$')
grid()
plot(f/1e9,s21_dB,'b-',linewidth=2, label='$S_{21}$')
legend()
ylabel('S-Parameter (dB)')
xlabel('Frequency (GHz)')

show()

kak = 1

