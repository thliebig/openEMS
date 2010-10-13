close all
clear
clc

%% setup the simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
length = 8000;
unit = 1e-3;
a = 1000;
width = a;
b = 500;
height = b;
mesh_res = [10 10 20];

%define mode
m = 1;
n = 0;

EPS0 = 8.85418781762e-12;
MUE0 = 1.256637062e-6;
C0 = 1/sqrt(EPS0*MUE0);

f0 = 1e9;

k = 2*pi*f0/C0;
kc = sqrt((m*pi/a/unit)^2 + (n*pi/b/unit)^2);
fc = C0*kc/2/pi;
beta = sqrt(k^2 - kc^2);

func_Ex = [num2str(n/b/unit) '*cos(' num2str(m*pi/a) '*x)*sin('  num2str(n*pi/b) '*y)'];
func_Ey = [num2str(m/a/unit) '*sin(' num2str(m*pi/a) '*x)*cos('  num2str(n*pi/b) '*y)'];

%% define and openEMS options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
openEMS_opts = '';
% openEMS_opts = [openEMS_opts ' --disable-dumps'];
% openEMS_opts = [openEMS_opts ' --debug-material'];
openEMS_opts = [openEMS_opts ' --engine=fastest'];

Sim_Path = 'tmp';
Sim_CSX = 'rect_wg.xml';

if (exist(Sim_Path,'dir'))
    rmdir(Sim_Path,'s');
end
mkdir(Sim_Path);

%% setup FDTD parameter & excitation function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FDTD = InitFDTD(50000,1e-5,'OverSampling',6);
FDTD = SetGaussExcite(FDTD,f0,f0/3);
BC = [0 0 0 0 2 2];
FDTD = SetBoundaryCond(FDTD,BC);

%% setup CSXCAD geometry & mesh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CSX = InitCSX();
mesh.x = 0 : mesh_res(1) : width;
mesh.y = 0 : mesh_res(2) : height;
mesh.z = 0 : mesh_res(3) : length;
CSX = DefineRectGrid(CSX, unit,mesh);

%% apply the excitation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
start=[0     0      mesh.z(3) ];
stop =[width height mesh.z(3) ];
CSX = AddExcitation(CSX,'excite',0,[1 1 0]);
weight{1} = func_Ex;
weight{2} = func_Ey;
weight{3} = 0;
CSX = SetExcitationWeight(CSX,'excite',weight);
CSX = AddBox(CSX,'excite',0 ,start,stop);
 
%% define dump boxes... %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CSX = AddDump(CSX,'Et','FileType',0,'SubSampling','4,4,1');
start = [mesh.x(1) , height/2 , mesh.z(1)];
stop = [mesh.x(end) ,  height/2 , mesh.z(end)];
CSX = AddBox(CSX,'Et',0 , start,stop);

% CSX = AddDump(CSX,'Ht','DumpType',1,'FileType',1,'SubSampling','4,4,4');
% CSX = AddBox(CSX,'Ht',0,start,stop);

%% Write openEMS compatoble xml-file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WriteOpenEMS([Sim_Path '/' Sim_CSX],FDTD,CSX);

RunOpenEMS(Sim_Path, Sim_CSX, openEMS_opts)

