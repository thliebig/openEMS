close all
clear
clc

%% setup the simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
abs_length = 250;
length = 1000;
width = 500;
height = 250;
MSL_width = 50;
MSL_height = 10;
mesh_res = [5 5 10];

EPS0 = 8.85418781762e-12;
MUE0 = 1.256637062e-6;

%% define openEMS options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
openEMS_opts = '';
% openEMS_opts = [openEMS_opts ' --disable-dumps'];
% openEMS_opts = [openEMS_opts ' --debug-material'];

Sim_Path = 'tmp';
Sim_CSX = 'msl.xml';

mkdir(Sim_Path);

%% setup FDTD parameter & excitation function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FDTD = InitFDTD(5e5,1e-6);
FDTD = SetGaussExcite(FDTD,0.5e9,0.5e9);
BC = [1 1 0 1 0 0];
FDTD = SetBoundaryCond(FDTD,BC);

%% setup CSXCAD geometry & mesh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CSX = InitCSX();
mesh.x = -width/2 : mesh_res(1) : width/2;
mesh.y = [linspace(0,MSL_height,11) MSL_height+1 MSL_height+3 MSL_height+mesh_res(2):mesh_res(2):height];
mesh.z = 0 : mesh_res(3) : length;
CSX = DefineRectGrid(CSX, 1e-3,mesh);

%%% MSL
CSX = AddMaterial(CSX,'copper');
CSX = SetMaterialProperty(CSX,'copper','Kappa',56e6);
start = [-0.5*MSL_width, MSL_height , 0];stop = [0.5*MSL_width, MSL_height+1 , length];
CSX = AddBox(CSX,'copper',0 ,start,stop);

start = [-0.5*MSL_width, 0 , 0];stop = [0.5*MSL_width, MSL_height , mesh_res(1)/2];
CSX = AddExcitation(CSX,'excite',0,[0 -1 0]);
CSX = AddBox(CSX,'excite',0 ,start,stop);

%% fake pml %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
finalKappa = 0.3/abs_length^4;
finalSigma = finalKappa*MUE0/EPS0;
CSX = AddMaterial(CSX,'pml');
CSX = SetMaterialProperty(CSX,'pml','Kappa',finalKappa);
CSX = SetMaterialProperty(CSX,'pml','Sigma',finalSigma);
CSX = SetMaterialWeight(CSX,'pml','Kappa',['pow(abs(z)-' num2str(length-abs_length) ',4)']);
CSX = SetMaterialWeight(CSX,'pml','Sigma',['pow(abs(z)-' num2str(length-abs_length) ',4)']);
start = [mesh.x(1) mesh.y(1) length-abs_length];
stop = [mesh.x(end) mesh.y(end) length];
CSX = AddBox(CSX,'pml',0 ,start,stop);
 
%% define dump boxes... %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CSX = AddDump(CSX,'Et_','DumpMode',2);
start = [mesh.x(1) , MSL_height/2 , mesh.z(1)];
stop = [mesh.x(end) , MSL_height/2 , mesh.z(end)];
CSX = AddBox(CSX,'Et_',0 , start,stop);

CSX = AddDump(CSX,'Ht_','DumpType',1,'DumpMode',2);
CSX = AddBox(CSX,'Ht_',0,start,stop);

%% define voltage calc boxes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%voltage calc
CSX = AddProbe(CSX,'ut1',0);
start = [ 0 MSL_height length/2 ];stop = [ 0 0 length/2 ];
CSX = AddBox(CSX,'ut1', 0 ,start,stop);

%current calc
CSX = AddProbe(CSX,'it1',1);
start = [ -MSL_width MSL_height/2 length/2 ];stop = [ MSL_width MSL_height*1.5 length/2 ];
CSX = AddBox(CSX,'it1', 0 ,start,stop);

%% Write openEMS compatoble xml-file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WriteOpenEMS([Sim_Path '/' Sim_CSX],FDTD,CSX);

%% cd to working dir and run openEMS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
savePath = pwd();
cd(Sim_Path); %cd to working dir
args = [Sim_CSX ' ' openEMS_opts];
invoke_openEMS(args);
cd(savePath);

