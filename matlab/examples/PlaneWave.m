close all
clear
clc

%% setup the simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
abs_length = 250;
length = 4000;
width = 1000;
height = 1000;
mesh_res = 25;

EPS0 = 8.85418781762e-12;
MUE0 = 1.256637062e-6;

%% define openEMS options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
openEMS_opts = '';
% openEMS_opts = [openEMS_opts ' --disable-dumps'];
% openEMS_opts = [openEMS_opts ' --debug-material'];
% openEMS_opts = [openEMS_opts ' --engine=multithreaded'];

Sim_Path = 'tmp';
Sim_CSX = 'plane_wave.xml';

mkdir(Sim_Path);

%% setup FDTD parameter & excitation function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FDTD = InitFDTD(5e5,1e-5,'OverSampling',10);
FDTD = SetGaussExcite(FDTD,0.5e9,0.5e9);
BC = [1 1 0 0 0 0];
FDTD = SetBoundaryCond(FDTD,BC);

%% setup CSXCAD geometry & mesh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CSX = InitCSX();
mesh.x = -width/2 : mesh_res : width/2;
mesh.y = -height/2 : mesh_res : height/2;
mesh.z = 0 : mesh_res : length;
CSX = DefineRectGrid(CSX, 1e-3,mesh);


%% fake pml %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
finalKappa = 0.3/abs_length^4;
finalSigma = finalKappa*MUE0/EPS0;
CSX = AddMaterial(CSX,'pml');
CSX = SetMaterialProperty(CSX,'pml','Kappa',finalKappa);
CSX = SetMaterialProperty(CSX,'pml','Sigma',finalSigma);
CSX = SetMaterialWeight(CSX,'pml','Kappa',['pow(abs(z)-' num2str(length-abs_length) ',4)']);
CSX = SetMaterialWeight(CSX,'pml','Sigma',['pow(abs(z)-' num2str(length-abs_length) ',4)']);
start=[-width/2 -height/2 length-abs_length];
stop=[width/2 height/2 length];
CSX = AddBox(CSX,'pml',0 ,start,stop);

%% apply the excitation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
start=[-width/2 -height/2 0];
stop=[width/2 height/2 0];
CSX = AddExcitation(CSX,'excite',0,[0 1 0]);
CSX = AddBox(CSX,'excite',0 ,start,stop);
 
%% define dump boxes... %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CSX = AddDump(CSX,'Et','FileType',1,'SubSampling','4,4,1');
start = [mesh.x(1) , mesh.y(1) , mesh.z(1)];
stop = [mesh.x(end) , mesh.y(end) , mesh.z(end)];
CSX = AddBox(CSX,'Et',0 , start,stop);

CSX = AddDump(CSX,'Ht','DumpType',1,'FileType',1,'SubSampling','4,4,1');
CSX = AddBox(CSX,'Ht',0,start,stop);

%% Write openEMS compatoble xml-file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WriteOpenEMS([Sim_Path '/' Sim_CSX],FDTD,CSX);

%% cd to working dir and run openEMS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
savePath = pwd();
cd(Sim_Path); %cd to working dir
args = [Sim_CSX ' ' openEMS_opts];
invoke_openEMS(args);
cd(savePath);

%% do the plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PlotArgs.slice = {mesh.x(round(end/2)) mesh.y(round(end/2)) mesh.z(round(end/2))};
PlotArgs.pauseTime=0.01;
PlotArgs.component=2;
PlotArgs.Limit = 'auto';

PlotHDF5FieldData('tmp/Et.h5',PlotArgs)
