close all;
clear all;
clc

feed_length=10;
wire_rad = sqrt(1.4/pi);
coil_rad = 10;
coil_length = 50;
coil_turns = 8;
coil_res = 10;
port_length = 10;
port_resist = 50;


openEMS_Path = [pwd() '/../../']
openEMS_opts = '';
% openEMS_opts = [openEMS_opts ' --disable-dumps'];
% openEMS_opts = [openEMS_opts ' --debug-material'];

Sim_Path = 'tmp';
Sim_CSX = 'helix.xml';

mkdir(Sim_Path);

%setup FDTD parameter
FDTD = InitFDTD(5e5,1e-6);
FDTD = SetGaussExcite(FDTD,0.5e9,0.5e9);
BC = [1 1 1 1 1 1];
FDTD = SetBoundaryCond(FDTD,BC);

%setup CSXCAD geometry
CSX = InitCSX();
mesh.x = [-35,-25,-20,-15,-12,-11,-10.5,-10,-9.5,-9,-8.5,-8,-7.5,-7,-6.5,-6,-5.5,-5,-4.5,-4,-3.5,-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10,10.5,11,12,13.5,15,17,18,19,19.5,20,20.5,21,22,23,25,27.5,30,35,45];
mesh.y = [-35,-25,-20,-15,-12,-11,-10.5,-10,-9.5,-9,-8.5,-8,-7.5,-7,-6.5,-6,-5.5,-5,-4.5,-4,-3.5,-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10,10.5,11,12,13,15,17.5,20,25,35];
mesh.z = [-25,-15,-10,-5,-2,-1,-0.5,0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10,10.5,11,11.5,12,12.5,13,13.5,14,14.5,15,15.5,16,16.5,17,17.5,18,18.5,19,19.5,20,20.5,21,21.5,22,22.5,23,23.5,24,24.5,25,25.5,26,26.5,27,27.5,28,28.5,29,29.5,30,30.5,31,31.5,32,32.5,33,33.5,34,34.5,35,35.5,36,36.5,37,37.5,38,38.5,39,39.5,40,40.5,41,41.5,42,42.5,43,43.5,44,44.5,45,45.5,46,46.5,47,47.5,48,48.5,49,49.5,50,50.5,51,52,53,55,57.5,60,65,75];
CSX = DefineRectGrid(CSX, 1e-3,mesh);

%create copper helix and feed lines...
CSX = AddMaterial(CSX,'copper');
CSX = SetMaterialProperty(CSX,'copper','Kappa',56e6);

%build helix-wire
dt = 1.0/coil_res;
height=0;
wire.Vertex = {};
count=0;
for n=0:coil_turns-1
    for m=0:coil_res
        count = count + 1;
        p(1,count) = coil_rad * cos(2*pi*dt*m);
        p(2,count) = coil_rad * sin(2*pi*dt*m);
        p(3,count) = height + coil_length/coil_turns * dt*m;
    end
    height = height + coil_length/coil_turns;
end
CSX = AddWire(CSX, 'copper', 0, p, wire_rad);

start = [coil_rad, 0 , 0];stop = [coil_rad+feed_length, 0 , 0];
CSX = AddCylinder(CSX,'copper',0 ,start,stop,wire_rad);

start(3)=coil_length;stop(3)=coil_length;
CSX = AddCylinder(CSX,'copper',0 ,start,stop,wire_rad);

start(3)=0;start(1)=coil_rad+feed_length;
CSX = AddCylinder(CSX,'copper',0 ,start,stop,wire_rad);

CSX = AddMaterial(CSX,'resist');
kappa = port_length/port_resist/wire_rad^2/pi/1e-3
CSX = SetMaterialProperty(CSX,'resist','Kappa',kappa);

start(3)=(coil_length-port_length)/2;stop(3)=(coil_length+port_length)/2;
CSX = AddCylinder(CSX,'resist',5 ,start,stop,wire_rad);

%excitation
CSX = AddExcitation(CSX,'excite',0,[0 0 1]);
CSX = AddCylinder(CSX,'excite', 0 ,start,stop,wire_rad);

%voltage calc
CSX = AddProbe(CSX,'ut1',0);
CSX = AddBox(CSX,'ut1', 0 ,start,stop);

%current calc
CSX = AddProbe(CSX,'it1',1);
start(3) = coil_length/2;stop(3) = coil_length/2;
start(1) = start(1)-2;start(2) = start(2)-2;
stop(1) = stop(1)+2;stop(2) = stop(2)+2;
CSX = AddBox(CSX,'it1', 0 ,start,stop);

%Write openEMS compatoble xml-file
WriteOpenEMS([Sim_Path '/' Sim_CSX],FDTD,CSX);

%cd to working dir and run openEMS
savePath = pwd();
cd(Sim_Path); %cd to working dir
command = [openEMS_Path 'openEMS ' Sim_CSX ' ' openEMS_opts];
disp(command);
system(command)
cd(savePath);

