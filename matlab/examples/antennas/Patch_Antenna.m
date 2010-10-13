close all
clear
clc

%% setup the simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% all length in mm
unit = 1e-3;

f0 = 2e9;
fc = 1e9;
freq = linspace(f0-fc,f0+fc,501);

physical_constants;

max_res = c0 / (f0+fc) / unit / 20;

% width in x-direction
% length in y-direction
% main radiation in z-direction
patch.width  = 32.86; %resonant length
patch.length = 41.37;

substrate.epsR   = 3.38;
substrate.width  = 120;
substrate.length = 120;
substrate.thickness = 1.524;
substrate.cells = 5;

feed.pos = -4.5;
feed.width = 0.5;
feed.R = 50; %feed resistance

%size of the simulation box
SimBox = [120 120 32];

%% define openEMS options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
openEMS_opts = '';
openEMS_opts = [openEMS_opts ' --engine=fastest'];
openEMS_opts = [openEMS_opts ' --numThreads=4'];

Sim_Path = 'tmp';
Sim_CSX = 'patch_ant.xml';

[status, message, messageid] = rmdir(Sim_Path,'s');
[status, message, messageid] = mkdir(Sim_Path);

%% setup FDTD parameter & excitation function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FDTD = InitFDTD(30000, 1e-5);
FDTD = SetGaussExcite(FDTD,f0,fc);
BC = [2 2 2 2 0 2]; %mur ABC
% BC = [3 3 3 3 0 3]; %use pml instead of mur
FDTD = SetBoundaryCond(FDTD,BC);

%% setup CSXCAD geometry & mesh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CSX = InitCSX();
mesh.x = [-SimBox(1)/2-8*max_res SimBox(1)/2+8*max_res -SimBox(1)/2 SimBox(1)/2 -substrate.width/2 substrate.width/2 feed.pos-feed.width/2 feed.pos+feed.width/2];
%add patch mesh with 2/3 - 1/3 rule
mesh.x = sort(unique([mesh.x -patch.width/2-max_res*0.66 -patch.width/2+max_res*0.33 patch.width/2+max_res*0.66 patch.width/2-max_res*0.33]));
mesh.x = SmoothMeshLines(mesh.x,max_res);
mesh.y = [-SimBox(2)/2-8*max_res SimBox(2)/2+8*max_res -SimBox(2)/2 SimBox(2)/2 -substrate.length/2 substrate.length/2 -feed.width/2 feed.width/2];
%add patch mesh with 2/3 - 1/3 rule
mesh.y = sort(unique([mesh.y -patch.length/2-max_res*0.66 -patch.length/2+max_res*0.33 patch.length/2+max_res*0.66 patch.length/2-max_res*0.33]));
mesh.y = SmoothMeshLines(mesh.y,max_res);
mesh.z = SmoothMeshLines([linspace(0,substrate.thickness,substrate.cells) SimBox(3) SimBox(3)+8*max_res],max_res);
CSX = DefineRectGrid(CSX, unit,mesh);

%% patch %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CSX = AddMetal(CSX,'patch');
start = [-patch.width/2 -patch.length/2 substrate.thickness];
stop  = [ patch.width/2  patch.length/2 substrate.thickness];
CSX = AddBox(CSX,'patch',10,start,stop);

%% substrate
CSX = AddMaterial(CSX,'substrate');
CSX = SetMaterialProperty(CSX,'substrate','Epsilon',substrate.epsR);
start = [-substrate.width/2 -substrate.length/2 0];
stop  = [ substrate.width/2  substrate.length/2 substrate.thickness];
CSX = AddBox(CSX,'substrate',0,start,stop);

%% apply the excitation & resist as a current source%%%%%%%%%%%%%%%%%%%%%%%
CSX = AddMaterial(CSX, 'resist');
kappa = substrate.thickness/feed.R/feed.width^2/unit;
CSX = SetMaterialProperty(CSX, 'resist', 'Kappa', kappa);
start=[feed.pos-feed.width/2 -feed.width/2 0];
stop =[feed.pos+feed.width/2  feed.width/2 substrate.thickness];
CSX = AddBox(CSX, 'resist', 15, start, stop);

CSX = AddExcitation(CSX, 'excite', 0, [0 0 1]);
CSX = AddBox(CSX, 'excite', 0, start, stop);

%% define voltage calc boxes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%voltage calc
start=[feed.pos 0 0];
stop =[feed.pos 0 substrate.thickness];
CSX = AddProbe(CSX,'ut1',0);
CSX = AddBox(CSX,'ut1', 0 ,stop,start);

%current calc
CSX = AddProbe(CSX,'it1',1);
start=[feed.pos-feed.width -feed.width substrate.thickness/2];
stop =[feed.pos+feed.width  feed.width substrate.thickness/2];
CSX = AddBox(CSX,'it1', 0 ,start,stop);

%% dump magnetic field over the patch antenna
CSX = AddDump(CSX,'Ht_','DumpType',1,'DumpMode',2);
start = [-patch.width -patch.length substrate.thickness+1];
stop  = [ patch.width  patch.length substrate.thickness+1];
CSX = AddBox(CSX,'Ht_',0 , start,stop);

%% Write openEMS compatoble xml-file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WriteOpenEMS([Sim_Path '/' Sim_CSX],FDTD,CSX);

%% run openEMS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RunOpenEMS(Sim_Path, Sim_CSX, openEMS_opts);

%% postproc & do the plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U = ReadUI('ut1','tmp/',freq);
I = ReadUI('it1','tmp/',freq);

close all

plot(U.TD{1}.t,U.TD{1}.val,'Linewidth',2);
grid on;

Zin = U.FD{1}.val./I.FD{1}.val;
figure()
plot(freq,real(Zin),'k-','Linewidth',2);
hold on;
grid on;
plot(freq,imag(Zin),'r--','Linewidth',2);

uf_inc = 0.5*(U.FD{1}.val + I.FD{1}.val * 50);
if_inc = 0.5*(I.FD{1}.val - U.FD{1}.val / 50);
uf_ref = U.FD{1}.val - uf_inc;
if_ref = I.FD{1}.val - if_inc;

s11 = uf_ref./uf_inc;
figure()
plot(freq,20*log10(abs(s11)),'k-','Linewidth',2);
grid on;

