close all;
clear all;
clc

abs_length = 500;
length = 5000;
unit = 1e-3;
rad  = 300;
mesh_res = [15 15 30];

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

func_Ex = [func_Er '*cos(a) - ' func_Ea '*sin(a)'];
func_Ey = [func_Er '*sin(a) + ' func_Ea '*cos(a)'];

openEMS_Path = [pwd() '/../../']
openEMS_opts = '';
% openEMS_opts = [openEMS_opts ' --disable-dumps'];
% openEMS_opts = [openEMS_opts ' --debug-material'];
% openEMS_opts = [openEMS_opts ' --debug-operator'];
% openEMS_opts = [openEMS_opts ' --engine=multithreaded'];

Sim_Path = 'tmp';
Sim_CSX = 'Circ_WG.xml';

mkdir(Sim_Path);

%setup FDTD parameter
FDTD = InitFDTD(1000,1e-6,'OverSampling',5);
T = 1/f0;
FDTD = SetCustomExcite(FDTD,f0,[ '(1-exp(-1*(t/' num2str(T) ')^2) ) * sin(2*pi*' num2str(f0) '*t)' ]);
BC = [1 1 1 1 1 1] * 0;
FDTD = SetBoundaryCond(FDTD,BC);

%setup CSXCAD geometry
CSX = InitCSX();
mesh.x = -mesh_res(1)/2-rad:mesh_res(1):rad+mesh_res(1)/2;
mesh.y = -mesh_res(2)/2-rad:mesh_res(2):rad+mesh_res(2)/2;
mesh.z = 0 : mesh_res(3) : length;
CSX = DefineRectGrid(CSX, 1e-3,mesh);

start = [0,0,0];
stop = [0,0,length];

%%fake pml
finalKappa = 0.3/abs_length^4;
finalSigma = finalKappa*MUE0/EPS0;
CSX = AddMaterial(CSX,'pml');
CSX = SetMaterialProperty(CSX,'pml','Kappa',finalKappa);
CSX = SetMaterialProperty(CSX,'pml','Sigma',finalSigma);
CSX = SetMaterialWeight(CSX,'pml','Kappa',['pow(abs(z)-' num2str(length-abs_length) ',4)']);
CSX = SetMaterialWeight(CSX,'pml','Sigma',['pow(abs(z)-' num2str(length-abs_length) ',4)']);
CSX = AddCylinder(CSX,'pml',10 ,[0 0 length-abs_length],stop,rad);

%%% fill everything with copper, priority 0
CSX = AddMaterial(CSX,'copper');
CSX = SetMaterialProperty(CSX,'copper','Kappa',56e6);
CSX = AddBox(CSX,'copper',0,[mesh.x(1) mesh.y(1) mesh.z(1)],[mesh.x(end) mesh.y(end) mesh.z(end)]);

%%% cut out an air cylinder as circular waveguide... priority 5
CSX = AddMaterial(CSX,'air');
CSX = SetMaterialProperty(CSX,'air','Epsilon',1);
CSX = AddCylinder(CSX,'air', 5 ,start,stop,rad);

CSX = AddExcitation(CSX,'excite',0,[1 1 0]);
weight{1} = func_Ex;
weight{2} = func_Ey;
weight{3} = 0;
CSX = SetExcitationWeight(CSX, 'excite', weight );
CSX = AddCylinder(CSX,'excite', 5 ,[0 0 -0.1],[0 0 0.1],rad);
 
%dump
CSX = AddDump(CSX,'Et','SubSampling','2,2,4','FileType',1,'DumpMode',2);
start = [mesh.x(1) , mesh.y(1) , mesh.z(1)];
stop = [mesh.x(end) , mesh.y(end) , mesh.z(end)];
CSX = AddBox(CSX,'Et',0 , start,stop);

% CSX = AddDump(CSX,'Ht','SubSampling','2,2,4','DumpType',1,'FileType',1,'DumpMode',2);
% CSX = AddBox(CSX,'Ht',0,start,stop);

% CSX = AddDump(CSX,'Excite_');
% start = [mesh.x(1) , mesh.y(1) , 0];
% stop = [mesh.x(end) , mesh.y(end) ,0];
% CSX = AddBox(CSX,'Excite_',0 , start,stop);
% 
% CSX = AddDump(CSX,'Exy');
% start = [mesh.x(1) , mesh.y(1) , length/2];
% stop = [mesh.x(end) , mesh.y(end) , length/2];
% CSX = AddBox(CSX,'Exy',0 , start,stop);

%voltage calc
CSX = AddProbe(CSX,'ut1',0);
start = [ -rad 0 0/2 ];stop = [ rad 0 0/2 ];
CSX = AddBox(CSX,'ut1', 0 ,start,stop);
% 
% %current calc
% CSX = AddProbe(CSX,'it1',1);
% mid = 0.5*(coax_rad_i+coax_rad_ai);
% start = [ -mid -mid length/2 ];stop = [ mid mid length/2 ];
% CSX = AddBox(CSX,'it1', 0 ,start,stop);

%Write openEMS compatoble xml-file
WriteOpenEMS([Sim_Path '/' Sim_CSX],FDTD,CSX);

%cd to working dir and run openEMS
savePath = pwd();
cd(Sim_Path); %cd to working dir
command = [openEMS_Path 'openEMS.sh ' Sim_CSX ' ' openEMS_opts];
disp(command);
system(command)
cd(savePath);

UI = ReadUI('ut1','tmp/');
plot(UI.TD{1}.t,UI.TD{1}.val);
grid on;

% plotting
if exist('tmp/Et.h5','file')
    PlotArgs.slice = {mesh.x(round(end/2)) mesh.y(round(end/2)) mesh.z(round(end/2))};
    PlotArgs.pauseTime=0.1;
    PlotArgs.component=0;
    PlotArgs.zlim='auto';

    PlotHDF5FieldData('tmp/Et.h5',PlotArgs)
end