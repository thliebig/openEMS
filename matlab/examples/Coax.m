close all;
clear all;
clc

abs_length = 250;
length = 1000;
coax_rad_i  = 100;
coax_rad_ai = 230;
coax_rad_aa = 240;
mesh_res = [5 5 5];

EPS0 = 8.85418781762e-12;
MUE0 = 1.256637062e-6;

openEMS_Path = [pwd() '/../../']
openEMS_opts = '';
% openEMS_opts = [openEMS_opts ' --disable-dumps'];
% openEMS_opts = [openEMS_opts ' --debug-material'];

Sim_Path = 'tmp';
Sim_CSX = 'coax.xml';

mkdir(Sim_Path);

%setup FDTD parameter
FDTD = InitFDTD(5e5,1e-6);
FDTD = SetGaussExcite(FDTD,0.5e9,0.5e9);
BC = [1 1 1 1 1 1] * 0;
FDTD = SetBoundaryCond(FDTD,BC);

%setup CSXCAD geometry
CSX = InitCSX();
mesh.x = -2.5*mesh_res(1)-coax_rad_aa : mesh_res(1) : coax_rad_aa+2.5*mesh_res(1);
mesh.y = mesh.x;
mesh.z = 0 : mesh_res(3) : length;
CSX = DefineRectGrid(CSX, 1e-3,mesh);

%create copper helix and feed lines...
CSX = AddMaterial(CSX,'copper');
CSX = SetMaterialProperty(CSX,'copper','Kappa',56e6);

%%%fake pml
finalKappa = 0.3/abs_length^4;
finalSigma = finalKappa*MUE0/EPS0;
CSX = AddMaterial(CSX,'pml');
CSX = SetMaterialProperty(CSX,'pml','Kappa',finalKappa);
CSX = SetMaterialProperty(CSX,'pml','Sigma',finalSigma);
CSX = SetMaterialWeight(CSX,'pml','Kappa',['pow(abs(z)-' num2str(length-abs_length) ',4)']);
CSX = SetMaterialWeight(CSX,'pml','Sigma',['pow(abs(z)-' num2str(length-abs_length) ',4)']);

%%% coax
start = [0, 0 , 0];stop = [0, 0 , length];
CSX = AddCylinder(CSX,'copper',0 ,start,stop,coax_rad_i);
CSX = AddCylindricalShell(CSX,'copper',0 ,start,stop,0.5*(coax_rad_aa+coax_rad_ai),(coax_rad_aa-coax_rad_ai));
start(3) = length-abs_length;
CSX = AddCylindricalShell(CSX,'pml',0 ,start,stop,0.5*(coax_rad_i+coax_rad_ai),(coax_rad_ai-coax_rad_i));
start(3) = 0; stop(3)=mesh_res(1)/2;
CSX = AddExcitation(CSX,'excite',0,[1 1 0]);
weight{1} = '(x)/(x*x+y*y)';
weight{2} = 'y/pow(rho,2)';
weight{3} = 0;
CSX = SetExcitationWeight(CSX, 'excite', weight );
CSX = AddCylindricalShell(CSX,'excite',0 ,start,stop,0.5*(coax_rad_i+coax_rad_ai),(coax_rad_ai-coax_rad_i));
 
%dump
CSX = AddDump(CSX,'Et_','DumpMode',2);
start = [mesh.x(1) , 0 , mesh.z(1)];
stop = [mesh.x(end) , 0 , mesh.z(end)];
CSX = AddBox(CSX,'Et_',0 , start,stop);

CSX = AddDump(CSX,'Ht_','DumpType',1,'DumpMode',2);
CSX = AddBox(CSX,'Ht_',0,start,stop);

%voltage calc
CSX = AddProbe(CSX,'ut1',0);
start = [ coax_rad_i 0 length/2 ];stop = [ coax_rad_ai 0 length/2 ];
CSX = AddBox(CSX,'ut1', 0 ,start,stop);

%current calc
CSX = AddProbe(CSX,'it1',1);
mid = 0.5*(coax_rad_i+coax_rad_ai);
start = [ -mid -mid length/2 ];stop = [ mid mid length/2 ];
CSX = AddBox(CSX,'it1', 0 ,start,stop);

%Write openEMS compatoble xml-file
WriteOpenEMS([Sim_Path '/' Sim_CSX],FDTD,CSX);

%cd to working dir and run openEMS
savePath = pwd();
cd(Sim_Path); %cd to working dir
command = [openEMS_Path 'openEMS.sh ' Sim_CSX ' ' openEMS_opts];
disp(command);
system(command)
cd(savePath);

