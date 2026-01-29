"""
 Coplanar waveguide with waveguide ports and
 local absorbers as terminations.

 (c) 2023-2025 Gadi Lahav <gadi@rfwithcare.com>

"""

### Import Libraries
import os, tempfile, shutil
from pylab import *
import scipy.io

from CSXCAD  import ContinuousStructure
from openEMS import openEMS
from openEMS.physical_constants import *
from CSXCAD.CSProperties import ABCtype

### General parameter setup
Sim_Path = os.path.join(tempfile.gettempdir(), 'Test_CPW_WG_Ports')

print(Sim_Path)

post_proc_only = False

# patch width (resonant length) in x-direction
Line_W = 1.15 #
CPW_gap = 0.3

#substrate setup
substrate_epsR   = 4.3
substrate_width  = 11
substrate_length = 50
substrate_thickness = 1

gap_cells = 2
substrate_cells = 3


# Port setup
port_w_fact = 7.5
port_h_fact = 3.5
portThick_mm = 0.5

cu_thick = 0.1

Airbox_Add = 12.5;
unit_res = 1e-3

# setup FDTD parameter & excitation function
f0 = 2e9 # center frequency
fc = 1e9 # 20 dB corner frequency


# size of the simulation box
SimBox = np.array([
            -substrate_width*0.5 - Airbox_Add,
            substrate_width*0.5 + Airbox_Add,
            -Airbox_Add,
            substrate_length + Airbox_Add,
            -substrate_thickness*(port_h_fact - 1.0) -Airbox_Add, 
            substrate_thickness*(1.0 + port_h_fact) + Airbox_Add])


print("Copying files ""CPW_E.csv"" and ""CPW_Hr.csv"" to {}".format(Sim_Path))
if not os.path.exists(Sim_Path):
    os.mkdir(Sim_Path)
shutil.copy("CPW_E.csv",Sim_Path)
shutil.copy("CPW_H.csv",Sim_Path)

### FDTD setup
## * Limit the simulation to 40k timesteps
## * Define a reduced end criteria of -40dB
FDTD = openEMS(NrTS=40000, EndCriteria=1e-4)
FDTD.SetGaussExcite( f0, fc )
FDTD.SetBoundaryCond( ['MUR', 'MUR', 'MUR', 'MUR', 'MUR', 'MUR'] )

CSX = ContinuousStructure()
FDTD.SetCSX(CSX)
mesh = CSX.GetGrid()
mesh.SetDeltaUnit(unit_res)
# mesh_res = ((C0/(f0+fc))/1e-3)/100
mesh_res = ((C0/(f0+fc))/1e-3)/50

### Generate propertiesprimitives and mesh-grid
#initialize the mesh with the "air-box" dimensions
mesh.AddLine('x', SimBox[0:2])
mesh.AddLine('y', SimBox[2:4])
mesh.AddLine('z', SimBox[4:6])

# create line and ground
line = CSX.AddMetal('cu_top')
start = [-Line_W/2, 0.0, substrate_thickness]
stop  = [ Line_W/2, substrate_length, substrate_thickness + cu_thick]
line.AddBox(priority = 20, start = start, stop = stop) # add a box-primitive to the metal property 'patch'
FDTD.AddEdges2Grid(dirs='xyz', properties = line)
mesh.AddLine('x',[start[0],stop[0]])
mesh.AddLine('y',[start[1],stop[1]])
mesh.AddLine('z',[start[2],stop[2]])

# PEC block for Z = 0
pec = CSX.AddMetal('pec')

start = [-Line_W*(1.0 + port_w_fact)*0.5, -portThick_mm, -substrate_thickness*(port_h_fact - 1.0)]
stop  = [ Line_W*(1.0 + port_w_fact)*0.5, 0.0          , substrate_thickness*(1.0 + port_h_fact) ]
# pec.AddBox(priority = 20, start = start, stop = stop)
mesh.AddLine('x',[start[0],stop[0]])
mesh.AddLine('y',[start[1],stop[1]])
mesh.AddLine('z',[start[2],stop[2]])

portStart = start
portStop = stop

# Extra pad for port not to take a lot of space
mesh.AddLine('y',np.array([0.125,0.25,0.375,1])*mesh_res)

start = [-Line_W*(1.0 + port_w_fact)*0.5, substrate_length               , -substrate_thickness*(port_h_fact - 1.0)]
stop  = [ Line_W*(1.0 + port_w_fact)*0.5, substrate_length + portThick_mm,  substrate_thickness*(1.0 + port_h_fact)]
# pec.AddBox(priority = 20, start = start, stop = stop)
mesh.AddLine('x',[start[0],stop[0]])
mesh.AddLine('y',[start[1],stop[1]])
mesh.AddLine('z',[start[2],stop[2]])

# Extra pad for port not to take a lot of space
mesh.AddLine('y',substrate_length - np.array([0.125,0.25,0.375,1])*mesh_res)

# create substrate
sub = CSX.AddMaterial('FR4', epsilon = substrate_epsR)
start = [-substrate_width/2, 0.0             , 0.0                ]
stop  = [ substrate_width/2, substrate_length, substrate_thickness]
sub.AddBox(priority = 2, start = start, stop = stop)
mesh.AddLine('x',[start[0],stop[0]])
mesh.AddLine('y',[start[1],stop[1]])
mesh.AddLine('z',[start[2],stop[2]])

# add extra cells to discretize the substrate thickness
mesh.AddLine('z', linspace(0, substrate_thickness, substrate_cells+1))

# create grounds
gnd = CSX.AddMetal('cu_gnd')
start = [-0.5*substrate_width     , 0.0             , substrate_thickness           ]
stop  = [-0.5*(Line_W + CPW_gap*2), substrate_length, substrate_thickness + cu_thick]
gnd.AddBox(priority=10, start=start, stop=stop) # add a box-primitive to the metal property 'patch'
mesh.AddLine('x',[start[0],stop[0]])
mesh.AddLine('y',[start[1],stop[1]])
mesh.AddLine('z',[start[2],stop[2]])


start = [0.5*(Line_W + CPW_gap*2), substrate_length, substrate_thickness + cu_thick]
stop  = [0.5*substrate_width     , 0.0             , substrate_thickness           ]
gnd.AddBox(priority=10, start=start, stop=stop) # add a box-primitive to the metal property 'patch'
mesh.AddLine('x',[start[0],stop[0]])
mesh.AddLine('y',[start[1],stop[1]])
mesh.AddLine('z',[start[2],stop[2]])

# Add extra mesh lines for CPW gaps
mesh.AddLine('x',linspace(-0.5*(Line_W + CPW_gap*2),-Line_W*0.5,gap_cells))
mesh.AddLine('x',linspace(Line_W*0.5,0.5*(Line_W + CPW_gap*2),gap_cells))

#################
# Generate mesh #
#################

mesh.SmoothMeshLines('all', mesh_res, 1.4)

# Find start\stop mesh lines
Yz = mesh.GetLines('y')
idxPort1 = (np.where(Yz == 0.0)[0] + 2).item(0)
idxPort2 = (np.where(Yz == substrate_length)[0] - 2).item(0)
idxAbs1 = idxPort1 - 1
idxAbs2 = idxPort2 + 1


# Calculate v_phase
beta = 62.5
v_phase = 2*pi*f0/beta


# Define port mode
start = [-Line_W*(1.0 + port_w_fact)*0.5, Yz.item(idxPort1 + 0), -substrate_thickness*(port_h_fact - 1.0)]
stop  = [ Line_W*(1.0 + port_w_fact)*0.5, Yz.item(idxPort1 + 1),  substrate_thickness*(1.0 + port_h_fact)] 
port1 = FDTD.AddWaveGuidePort( 1, start, stop, 'y', E_file = "CPW_E.csv", H_file = "CPW_H.csv", kc = 0.0, excite = 1, excite_type = 0)

# Place absorber for port #1
stop[1]  = Yz.item(idxAbs1)
start[1] = stop[1] # Place absorber at end of port
abs1 = CSX.AddAbsorbingBC('abs1',NormalSignPositive = True, AbsorbingBoundaryType = ABCtype.MUR_1ST)
abs1.AddBox(start, stop, priority = 20)
abs1.SetPhaseVelocity(v_phase);

start = [-Line_W*(1.0 + port_w_fact)*0.5, Yz.item(idxPort2 - 0), -substrate_thickness*(port_h_fact - 1.0)]
stop  = [ Line_W*(1.0 + port_w_fact)*0.5, Yz.item(idxPort2 - 1),  substrate_thickness*(1.0 + port_h_fact)]
port2 = FDTD.AddWaveGuidePort( 2, start, stop, 'y', E_file = "CPW_E.csv", H_file = "CPW_H.csv", kc = 0.0, excite = 0, excite_type = 0)

# Place absorber for port #2
start[1] = Yz.item(idxAbs2)
stop[1]  = start[1]
abs2 = CSX.AddAbsorbingBC('abs2',NormalSignPositive = False, AbsorbingBoundaryType = ABCtype.MUR_1ST)
abs2.AddBox(start, stop, priority=20)
abs2.SetPhaseVelocity(v_phase);

# # Define dump box...
# Et = CSX.AddDump('Et', file_type=0, dump_type=0, dump_mode=1)
# start = [float(SimBox[0]), float(SimBox[2]), float(SimBox[4])];
# stop  = [float(SimBox[1]), float(SimBox[3]), float(SimBox[5])];
# Et.AddBox(start, stop);

### Run the simulation
if 1:  # debugging only
    CSX_file = os.path.join(Sim_Path, 'cpw_2_WG_ports.xml')
    if not os.path.exists(Sim_Path):
        os.mkdir(Sim_Path)
    CSX.Write2XML(CSX_file)
    from CSXCAD import AppCSXCAD_BIN
    os.system(AppCSXCAD_BIN + ' "{}"'.format(CSX_file))

if not post_proc_only:
    FDTD.Run(Sim_Path, verbose=0, cleanup=False)

Zl = 45
Zw = 254

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
