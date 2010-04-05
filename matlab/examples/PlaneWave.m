close all;
clear all;
clc

abs_length = 250;
length = 4000;
width = 1000;
height = 1000;
mesh_res = 25;

EPS0 = 8.85418781762e-12;
MUE0 = 1.256637062e-6;

openEMS_Path = [pwd() '/../../']
openEMS_opts = '';
% openEMS_opts = [openEMS_opts ' --disable-dumps'];
% openEMS_opts = [openEMS_opts ' --debug-material'];

Sim_Path = 'tmp';
Sim_CSX = 'plane_wave.xml';

mkdir(Sim_Path);

%setup FDTD parameter
FDTD = InitFDTD(5e5,1e-6);
FDTD = SetGaussExcite(FDTD,0.5e9,0.5e9);
BC = [1 1 0 0 0 0];
FDTD = SetBoundaryCond(FDTD,BC);

%setup CSXCAD geometry
CSX = InitCSX();
mesh.x = -width/2 : mesh_res : width/2;
mesh.y = -height/2 : mesh_res : height/2;
mesh.z = 0 : mesh_res : length;
CSX = DefineRectGrid(CSX, 1e-3,mesh);


%%%fake pml
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

start=[-width/2 -height/2 0];
stop=[width/2 height/2 0];
CSX = AddExcitation(CSX,'excite',0,[0 1 0]);
CSX = AddBox(CSX,'excite',0 ,start,stop);
 
%dump
CSX = AddDump(CSX,'Et',0,0,1);
start = [mesh.x(1) , 0 , mesh.z(1)];
stop = [mesh.x(end) , 0 , mesh.z(end)];
CSX = AddBox(CSX,'Et',0 , start,stop);

CSX = AddDump(CSX,'Ht',1,0,1);
CSX = AddBox(CSX,'Ht',0,start,stop);

%Write openEMS compatoble xml-file
WriteOpenEMS([Sim_Path '/' Sim_CSX],FDTD,CSX);

%cd to working dir and run openEMS
savePath = pwd();
cd(Sim_Path); %cd to working dir
command = [openEMS_Path 'openEMS.sh ' Sim_CSX ' ' openEMS_opts];
disp(command);
system(command)
cd(savePath);

