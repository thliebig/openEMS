close all;
clear all;
clc

abs_length = 500;
length = 5000;
unit = 1e-3;
rad  = 300;
mesh_max = 15;
N_alpha = ceil(rad * 2*pi / mesh_max);
mesh_res = [mesh_max 2*pi/N_alpha mesh_max];


EPS0 = 8.85418781762e-12;
MUE0 = 1.256637062e-6;
C0 = 1/sqrt(EPS0*MUE0);

f0 = 400e6;

p11 = 1.841;
kc = p11 / rad /unit;
k = 2*pi*f0/C0;
fc = C0*kc/2/pi
beta = sqrt(k^2 - kc^2);

kc = kc*unit;
func_Er = [ num2str(-1/kc^2) '/rho*cos(a)*j1('  num2str(kc) '*rho)'];
func_Ea = [ num2str(1/kc) '*sin(a)*0.5*(j0('  num2str(kc) '*rho)-jn(2,'  num2str(kc) '*rho))'];

openEMS_Path = [pwd() '/../../']
openEMS_opts = '';
% openEMS_opts = [openEMS_opts ' --disable-dumps'];
% openEMS_opts = [openEMS_opts ' --debug-material'];
% openEMS_opts = [openEMS_opts ' --debug-operator'];
% openEMS_opts = [openEMS_opts ' --engine=multithreaded'];

Sim_Path = 'tmp';
Sim_CSX = 'Circ_WG_CC.xml';

mkdir(Sim_Path);

%setup FDTD parameter
FDTD = InitCylindricalFDTD(1e5,1e-5,'OverSampling',10);
T = 1/f0;
FDTD = SetCustomExcite(FDTD,f0,[ '(1-exp(-1*(t/' num2str(T) ')^2) ) * sin(2*pi*' num2str(f0) '*t)' ]);
BC = [0 0 0 0 0 0];
FDTD = SetBoundaryCond(FDTD,BC);

%setup CSXCAD geometry
CSX = InitCSX();
mesh.x = [0 2*mesh_res(1):mesh_res(1):rad];
mesh.y = linspace(-pi,pi,N_alpha);
mesh.z = 0 : mesh_res(3) : length;
CSX = DefineRectGrid(CSX, 1e-3,mesh);

start = [0 mesh.y(1) length-abs_length];
stop = [rad mesh.y(end) length];

%%fake pml
finalKappa = 0.3/abs_length^4;
finalSigma = finalKappa*MUE0/EPS0;
CSX = AddMaterial(CSX,'pml');
CSX = SetMaterialProperty(CSX,'pml','Kappa',finalKappa);
CSX = SetMaterialProperty(CSX,'pml','Sigma',finalSigma);
CSX = SetMaterialWeight(CSX,'pml','Kappa',['pow(abs(z)-' num2str(length-abs_length) ',4)']);
CSX = SetMaterialWeight(CSX,'pml','Sigma',['pow(abs(z)-' num2str(length-abs_length) ',4)']);
CSX = AddBox(CSX,'pml',0 ,start,stop);

CSX = AddExcitation(CSX,'excite',0,[1 1 0]);
weight{1} = func_Er;
weight{2} = func_Ea;
weight{3} = 0;
CSX = SetExcitationWeight(CSX, 'excite', weight );
start(3)=-5;
stop(3)=5;
CSX = AddBox(CSX,'excite', 5 ,start,stop);
 
%dump
CSX = AddDump(CSX,'Et','FileType',0,'DumpMode',0);
start = [mesh.x(1) ,0 , mesh.z(1)];
stop = [mesh.x(end) , 0 , mesh.z(end)];
CSX = AddBox(CSX,'Et',0 , start,stop);

CSX = AddDump(CSX,'Ht','DumpType',1,'FileType',0,'DumpMode',0);
start = [mesh.x(1) ,0 , mesh.z(1)];
stop = [mesh.x(end) , 0 , mesh.z(end)];
CSX = AddBox(CSX,'Ht',0 , start,stop);

%voltage calc
CSX = AddProbe(CSX,'ut_exc',0);
start = [ 0 0 0 ];stop = [ rad 0 0 ];
CSX = AddBox(CSX,'ut_exc', 0 ,start,stop);

CSX = AddProbe(CSX,'ut_1',0);
start = [ 0 0 length/2 ];stop = [ rad 0 length/2 ];
CSX = AddBox(CSX,'ut_1', 0 ,start,stop);

%Write openEMS compatoble xml-file
WriteOpenEMS([Sim_Path '/' Sim_CSX],FDTD,CSX);

%cd to working dir and run openEMS
savePath = pwd();
cd(Sim_Path); %cd to working dir
command = [openEMS_Path 'openEMS.sh ' Sim_CSX ' ' openEMS_opts];
disp(command);
system(command)
cd(savePath);

UI = ReadUI('ut_1','tmp/');
plot(UI.TD{1}.t,UI.TD{1}.val);
grid on;

% plotting
% if exist('tmp/Et.h5','file')
%     PlotArgs.slice = {mesh.x(round(end/2)) mesh.y(round(end/2)) mesh.z(round(end/2))};
%     PlotArgs.pauseTime=0.1;
%     PlotArgs.component=0;
%     PlotArgs.zlim='auto';
% 
%     PlotHDF5FieldData('tmp/Et.h5',PlotArgs)
% end