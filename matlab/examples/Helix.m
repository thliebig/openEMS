close all
clear
clc

%% setup the simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
feed_length=10;
wire_rad = sqrt(1.4/pi);
mesh_size = wire_rad;
coil_rad = 10;
coil_length = 50;
coil_turns = 8;
coil_res = 10;
port_length = mesh_size; %coil_length/2;
port_resist = 1000;

f_max = 100e6;
f_excite = 300e6;

%% define openEMS options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
openEMS_opts = '';
% openEMS_opts = [openEMS_opts ' --debug-material'];
% openEMS_opts = [openEMS_opts ' --debug-boxes'];
% openEMS_opts = [openEMS_opts ' --debug-operator'];

openEMS_opts = [openEMS_opts ' --disable-dumps --engine=fastest'];
% openEMS_opts = [openEMS_opts ' --engine=sse-compressed'];

Sim_Path = 'tmp';
Sim_CSX = 'helix.xml';

[status, message, messageid] = rmdir(Sim_Path,'s');
[status,message,messageid] = mkdir(Sim_Path);

%% setup FDTD parameter & excitation function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FDTD = InitFDTD(30000,1e-6);
FDTD = SetGaussExcite(FDTD,f_excite/2,f_excite/2);
BC = [1 1 1 1 1 1];
FDTD = SetBoundaryCond(FDTD,BC);

%% setup CSXCAD geometry & mesh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
add_Lines = mesh_size * 1.5.^(1:10);
add_Lines = add_Lines(find(add_Lines<(3e8/f_excite)/10*1e3));

CSX = InitCSX();
mesh.x = -coil_rad-mesh_size : mesh_size : coil_rad+mesh_size+feed_length;
mesh.x = [mesh.x(1)-add_Lines mesh.x mesh.x(end)+add_Lines ];
mesh.y = -coil_rad-mesh_size : mesh_size : coil_rad+mesh_size;
mesh.y = [mesh.y(1)-add_Lines mesh.y mesh.y(end)+add_Lines ];
mesh.z = -mesh_size : mesh_size : coil_length+mesh_size;
mesh.z = [mesh.z(1)-add_Lines mesh.z mesh.z(end)+add_Lines ];
CSX = DefineRectGrid(CSX, 1e-3,mesh);

%% build/define helix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CSX = AddMaterial(CSX,'copper');
CSX = SetMaterialProperty(CSX,'copper','Kappa',56e6);

dt = 1.0/coil_res;
height=0;
wire.Vertex = {};
p(1,1) = coil_rad + feed_length;
p(2,1) = 0;
p(3,1) = 0.5*(coil_length-port_length);
p(1,2) = coil_rad + feed_length;
p(2,2) = 0;
p(3,2) = 0;
count=2;
for n=0:coil_turns-1
    for m=0:coil_res
        count = count + 1;
        p(1,count) = coil_rad * cos(2*pi*dt*m);
        p(2,count) = coil_rad * sin(2*pi*dt*m);
        p(3,count) = height + coil_length/coil_turns * dt*m;
    end
    height = height + coil_length/coil_turns;
end
p(1,count+1) = coil_rad + feed_length;
p(2,count+1) = 0;
p(3,count+1) = coil_length;
p(1,count+2) = coil_rad + feed_length;
p(2,count+2) = 0;
p(3,count+2) = 0.5*(coil_length+port_length);
CSX = AddWire(CSX, 'copper', 0, p, wire_rad);

%% apply the excitation & resist as a current source%%%%%%%%%%%%%%%%%%%%%%%
CSX = AddMaterial(CSX,'resist');
kappa = port_length/port_resist/wire_rad^2/pi/1e-3;
CSX = SetMaterialProperty(CSX,'resist','Kappa',kappa);

start=[coil_rad+feed_length 0 (coil_length-port_length)/2];
stop=[coil_rad+feed_length 0 (coil_length+port_length)/2];
%start(3)=(coil_length-port_length)/2;stop(3)=(coil_length+port_length)/2;
CSX = AddCylinder(CSX,'resist',5 ,start,stop,wire_rad);

CSX = AddExcitation(CSX,'excite',0,[0 0 1]);
CSX = AddCylinder(CSX,'excite', 0 ,start,stop,wire_rad);

%% define voltage calc boxes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%voltage calc
CSX = AddProbe(CSX,'ut1',0);
CSX = AddBox(CSX,'ut1', 0 ,stop,start);

%current calc
CSX = AddProbe(CSX,'it1',1);
start(3) = coil_length/2+mesh_size;stop(3) = coil_length/2+mesh_size;
start(1) = start(1)-2;start(2) = start(2)-2;
stop(1) = stop(1)+2;stop(2) = stop(2)+2;
CSX = AddBox(CSX,'it1', 0 ,start,stop);

%% define dump boxes... %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CSX = AddDump(CSX,'Et_');
start = [mesh.x(1) , 0 , mesh.z(1)];
stop = [mesh.x(end) , 0 , mesh.z(end)];
CSX = AddBox(CSX,'Et_',0 , start,stop);

CSX = AddDump(CSX,'Ht_','DumpType',1);
start = [mesh.x(1) , 0 , mesh.z(1)];
stop = [mesh.x(end) , 0 , mesh.z(end)];
CSX = AddBox(CSX,'Ht_',0 , start,stop);

%% Write openEMS compatoble xml-file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WriteOpenEMS([Sim_Path '/' Sim_CSX],FDTD,CSX);

%% run openEMS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RunOpenEMS(Sim_Path, Sim_CSX, openEMS_opts);

%% postproc & do the plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U = ReadUI('ut1','tmp/');
I = ReadUI('it1','tmp/');

Z = U.FD{1}.val./I.FD{1}.val;
f = U.FD{1}.f;
L = imag(Z)./(f*2*pi);
R = real(Z);
ind = find(f<f_max);

subplot(2,1,1);
plot(f(ind)*1e-6,L(ind)*1e9,'Linewidth',2);
xlabel('frequency (MHz)');
ylabel('coil inductance (nH)');
grid on;
subplot(2,1,2);
plot(f(ind)*1e-6,R(ind),'Linewidth',2);
hold on
plot(f(ind)*1e-6,imag(Z(ind)),'r','Linewidth',2);
xlabel('frequency (MHz)');
ylabel('resistance (Ohm)');
grid on;
legend( {'real','imaginary'}, 'location', 'northwest' )

figure
plot(U.TD{1}.t/1e-6,U.TD{1}.val,'Linewidth',2);
xlabel('time (us)');
ylabel('amplitude (V)');
grid on;
