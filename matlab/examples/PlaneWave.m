close all
clear
clc

%% setup the simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
length = 5000;
width = 300;
height = 200;
mesh_res = 15;
abs_length = mesh_res*10;

EPS0 = 8.85418781762e-12;
MUE0 = 1.256637062e-6;

%% define openEMS options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
openEMS_opts = '';
% openEMS_opts = [openEMS_opts ' --disable-dumps'];
% openEMS_opts = [openEMS_opts ' --debug-material'];
openEMS_opts = [openEMS_opts ' --engine=fastest'];

Sim_Path = 'tmp';
Sim_CSX = 'plane_wave.xml';

if (exist(Sim_Path,'dir'))
    rmdir(Sim_Path,'s');
end
mkdir(Sim_Path);

%% setup FDTD parameter & excitation function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FDTD = InitFDTD(5000,1e-5,'OverSampling',10);
FDTD = SetGaussExcite(FDTD,0.5e9,0.5e9);
BC = [1 1 0 0 2 2];
FDTD = SetBoundaryCond(FDTD,BC);

%% setup CSXCAD geometry & mesh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CSX = InitCSX();
mesh.x = -width/2 : mesh_res : width/2;
mesh.y = -height/2 : mesh_res : height/2;
mesh.z = 0 : mesh_res : length;
CSX = DefineRectGrid(CSX, 1e-3,mesh);

%% apply the excitation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
start=[-width/2 -height/2 mesh.z(3)];
stop =[ width/2  height/2 mesh.z(3)];
CSX = AddExcitation(CSX,'excite',0,[0 1 0]);
CSX = AddBox(CSX,'excite',0 ,start,stop);
 
%% define dump boxes... %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CSX = AddDump(CSX,'Et','FileType',1,'SubSampling','4,4,4');
start = [mesh.x(1) , mesh.y(1) , mesh.z(1)];
stop = [mesh.x(end) , mesh.y(end) , mesh.z(end)];
CSX = AddBox(CSX,'Et',0 , start,stop);

CSX = AddDump(CSX,'Ht','DumpType',1,'FileType',1,'SubSampling','4,4,4','DumpMode',2);
CSX = AddBox(CSX,'Ht',0,start,stop);

%% Write openEMS compatoble xml-file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WriteOpenEMS([Sim_Path '/' Sim_CSX],FDTD,CSX);

%% run openEMS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RunOpenEMS(Sim_Path, Sim_CSX, openEMS_opts);

%% do the plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PlotArgs.slice = {mesh.x(round(end/2)) mesh.y(round(end/2)) mesh.z(round(end/2))};
PlotArgs.pauseTime=0.01;
PlotArgs.component=1;
PlotArgs.Limit = 'auto';

PlotHDF5FieldData('tmp/Ht.h5',PlotArgs)
