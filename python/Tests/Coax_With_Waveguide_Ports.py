### Import Libraries
import os, tempfile
from pylab import *

from CSXCAD  import ContinuousStructure
from openEMS import openEMS
from openEMS.physical_constants import *
from CSXCAD.CSProperties import BCtype

import scipy.io

### General parameter setup
Sim_Path = os.path.join(tempfile.gettempdir(), 'Test_Coax_WG_Ports')

print(Sim_Path)

post_proc_only = False

# patch width (resonant length) in x-direction
microstrip_W = 1.875 #

#substrate setup
coax_D = 2
coax_shield_thick = 0.15
coax_wire_D = 0.5
coax_L = 25

teflon_epsR = 2.5

mesh_res = 0.5

# Port setup
port_w_fact = 6
port_h_fact = 5.5
portThick_mm = 0.5

cu_thick = 0.1

Airbox_Add = 1;

unit_res = 1e-3


# size of the simulation box
SimBox = np.array([
            -(coax_D*0.5 + coax_shield_thick + Airbox_Add),
            (coax_D*0.5 + coax_shield_thick + Airbox_Add),
            -(coax_D*0.5 + coax_shield_thick + Airbox_Add),
            (coax_D*0.5 + coax_shield_thick + Airbox_Add),
            -Airbox_Add, 
            coax_L + Airbox_Add])

# setup FDTD parameter & excitation function
f0 = 2.5e9 # center frequency
fc = 1e9 # 20 dB corner frequency

### FDTD setup
## * Limit the simulation to 30k timesteps
## * Define a reduced end criteria of -40dB
FDTD = openEMS(NrTS=30000, EndCriteria=1e-4)
FDTD.SetGaussExcite( f0, fc )
FDTD.SetBoundaryCond( ['MUR', 'MUR', 'MUR', 'MUR', 'MUR', 'MUR'] )

CSX = ContinuousStructure()
FDTD.SetCSX(CSX)

mesh = CSX.GetGrid()
mesh.SetDeltaUnit(unit_res)
mesh_res = ((C0/(f0+fc))/1e-3)/100
### Generate propertiesprimitives and mesh-grid
#initialize the mesh with the "air-box" dimensions
mesh.AddLine('x', SimBox[0:2])
mesh.AddLine('y', SimBox[2:4])
mesh.AddLine('z', SimBox[4:6])

# create center wire
line = CSX.AddMetal('Wire_Inner')
start = [0.0, 0.0, 0.0]
stop  = [0.0, 0.0, coax_L]
line.AddCylinder(priority=10,start=start,stop=stop,radius=coax_wire_D*0.5)
mesh.AddLine('x',np.linspace(start[0]-coax_wire_D*0.5,stop[0]+coax_wire_D*0.5,6).tolist())
mesh.AddLine('y',np.linspace(start[1]-coax_wire_D*0.5,stop[1]+coax_wire_D*0.5,6).tolist())
mesh.AddLine('z',[start[2],stop[2]])

# create Outer Shield
shield = CSX.AddMetal('Shield_Outer')
start = [0.0, 0.0, 0.0]
stop  = [0.0, 0.0, coax_L]
shield.AddCylindricalShell(priority=10,start=start,stop=stop,radius=(coax_D + coax_shield_thick)*0.5,shell_width=coax_shield_thick)
rad = (coax_D + coax_shield_thick)*0.5
hthick = coax_shield_thick*0.5
mesh.AddLine('x',
                np.linspace(start[0] - (rad + hthick),start[0] - (rad - hthick),4).tolist() + 
                np.linspace(stop[0]  + (rad - hthick),stop[0]  + (rad + hthick),4).tolist())
mesh.AddLine('y',
                np.linspace(start[1] - (rad + hthick),start[1] - (rad - hthick),4).tolist() + 
                np.linspace(stop[1]  + (rad - hthick),stop[1]  + (rad + hthick),4).tolist())
mesh.AddLine('z',[start[2],stop[2]])

# Create teflon fill
teflon = CSX.AddMaterial('PTFE',epsilon=teflon_epsR)
start = [0.0, 0.0, 0.0]
stop  = [0.0, 0.0, coax_L]
teflon.AddCylindricalShell(priority=8,start=start,stop=stop,radius=(coax_wire_D + coax_D)*0.25,shell_width=(coax_D - coax_wire_D)*0.5)
rad = (coax_wire_D + coax_D)*0.25
hthick = (coax_D - coax_wire_D)*0.25
mesh.AddLine('x',
                np.linspace(start[0] - (rad + hthick),start[0] - (rad - hthick),12).tolist() + 
                np.linspace(stop[0]  + (rad - hthick),stop[0]  + (rad + hthick),12).tolist())
mesh.AddLine('y',
                np.linspace(start[1] - (rad + hthick),start[1] - (rad - hthick),12).tolist() + 
                np.linspace(stop[1]  + (rad - hthick),stop[1]  + (rad + hthick),12).tolist())
mesh.AddLine('z',[start[2],stop[2]])

# Port PEC blocks
port_blocks = CSX.AddMetal('PEC')

start = [-(coax_D*0.5 + coax_shield_thick), -(coax_D*0.5 + coax_shield_thick), -portThick_mm]
stop = [(coax_D*0.5 + coax_shield_thick), (coax_D*0.5 + coax_shield_thick), 0.0]
port_blocks.AddBox(priority=12, start=start, stop=stop) 
mesh.AddLine('z',[-portThick_mm])

start = [-(coax_D*0.5 + coax_shield_thick), -(coax_D*0.5 + coax_shield_thick), coax_L]
stop = [(coax_D*0.5 + coax_shield_thick), (coax_D*0.5 + coax_shield_thick), coax_L + portThick_mm]
port_blocks.AddBox(priority=12, start=start, stop=stop)
mesh.AddLine('z',[coax_L + portThick_mm])

# Add dense mesh lines close to ports
mesh.AddLine('z',np.array([0.25,0.5,0.8,1])*mesh_res)
mesh.AddLine('z',coax_L - np.array([0.25,0.5,0.8,1])*mesh_res)

mesh.SmoothMeshLines('all', mesh_res, 1.25)

Zz = mesh.GetLines('z')
# Find start\stop mesh lines
idxPort1 = (np.where(Zz == 0.0)[0] + 2).item(0)
idxPort2 = (np.where(Zz == coax_L)[0] - 2).item(0)
idxAbs1 = idxPort1 - 1
idxAbs2 = idxPort2 + 1


matFile = scipy.io.loadmat("Coax_Mode_1.mat")
Em1 = matFile["E"]
Hm1 = matFile["H"]

# Define port mode 

start = [-coax_D*0.5-coax_shield_thick, -coax_D*0.5-coax_shield_thick, Zz.item(idxPort1 + 0)]
stop  = [coax_D*0.5+coax_shield_thick, coax_D*0.5+coax_shield_thick, Zz.item(idxPort1 + 1)]
port1 = FDTD.AddWaveGuidePort( 1, start, stop, 'z', Em1, Hm1, 0.0, excite=1, excite_type = 0)

# Place absorber for port #1
stop[2] = Zz.item(idxAbs1)
start[2] = stop[2] # Place absorber at end of port
abs1 = CSX.AddAbsorbingBC('abs1',NormDir = 3, BCtype = BCtype.MUR_1ST_1PV_SA)
abs1.AddBox(start, stop, priority=20)

# Place absorber for port #2
start[2] = Zz.item(idxAbs2)
stop[2] = start[2]
abs2 = CSX.AddAbsorbingBC('abs2',NormDir = -3, BCtype = BCtype.MUR_1ST_1PV_SA)
abs2.AddBox(start, stop, priority=20)

matFile = scipy.io.loadmat("Coax_Mode_1.mat")
Em2 = matFile["E"]
Hm2 = matFile["H"]
kz = matFile["Kz"]

start = [-coax_D*0.5-coax_shield_thick, -coax_D*0.5-coax_shield_thick, Zz.item(idxPort2 - 0)]
stop  = [coax_D*0.5+coax_shield_thick, coax_D*0.5+coax_shield_thick, Zz.item(idxPort2 - 1)]
port2 = FDTD.AddWaveGuidePort( 2, start, stop, 'z', Em2, Hm2, 0.0, excite=0, excite_type = 0)

v_phase = 2*np.pi*f0/kz
v_phase = v_phase[0]
abs1.SetPhaseVelocity(v_phase);
abs2.SetPhaseVelocity(v_phase);

# # # Define dump box...
# # Et = CSX.AddDump('Et', file_type=0, dump_type=0, dump_mode=1)
# # start = [float(SimBox[0]), float(SimBox[2]), float(SimBox[4])];
# # stop  = [float(SimBox[1]), float(SimBox[3]), float(SimBox[5])];
# # Et.AddBox(start, stop);

### Run the simulation
if 0: # debugging only
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
port1.CalcPort(Sim_Path, f,ref_impedance = Zw, ZL = Zl)
port2.CalcPort(Sim_Path, f,ref_impedance = Zw, ZL = Zl)
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





