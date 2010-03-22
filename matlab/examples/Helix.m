close all;
clear all;
clc

openEMS_Path = [pwd() '/../../']
openEMS_opts = '';
openEMS_opts = [openEMS_opts ' --disable-dumps'];
openEMS_opts = [openEMS_opts ' --debug-material'];

Sim_Path = 'tmp';
Sim_CSX = 'helix.xml';

mkdir(Sim_Path);

%setup FDTD parameter
FDTD = InitFDTD(1000,1e-6);
FDTD = SetGaussExcite(FDTD,0.5e9,0.5e9);
BC = [1 1 1 1 1 1];
FDTD = SetBoundaryCond(FDTD,BC);

%setup CSXCAD geometry
CSX = [];

%Write openEMS compatoble xml-file
WriteOpenEMS([Sim_Path '/' Sim_CSX],FDTD,CSX);

%cd to working dir and run openEMS
savePath = pwd();
cd(Sim_Path); %cd to working dir
command = [openEMS_Path 'openEMS ' openEMS_opts];
disp(command);
system(command)
cd(savePath);

