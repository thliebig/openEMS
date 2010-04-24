close all;
clear all;
clc

%% setup the simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EPS0 = 8.85418781762e-12;
MUE0 = 1.256637062e-6;
C0 = 1/sqrt(EPS0*MUE0);
Z0 = sqrt(MUE0/EPS0);

f0 = 0.5e9;
epsR = 1;

abs_length = 500;
length = 3000;
port_dist = 1500;
rad_i  = 100;
rad_a = 230;
max_mesh = 10;
max_alpha = max_mesh;
N_alpha = ceil(rad_a * 2*pi / max_alpha);
mesh_res = [max_mesh 2*pi/N_alpha max_mesh];

%% define file pathes and openEMS options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
openEMS_Path = [pwd() '/../../'];
openEMS_opts = '';
openEMS_opts = [openEMS_opts ' --disable-dumps'];
% openEMS_opts = [openEMS_opts ' --debug-material'];

Sim_Path = 'tmp';
Sim_CSX = 'coax.xml';

mkdir(Sim_Path);

%% setup FDTD parameter & excitation function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FDTD = InitCylindricalFDTD(1e5,1e-5,'OverSampling',10);
FDTD = SetGaussExcite(FDTD,f0,f0);
BC = [0 0 1 1 0 0];
FDTD = SetBoundaryCond(FDTD,BC);

%% setup CSXCAD geometry & mesh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CSX = InitCSX();
mesh.x = rad_i : mesh_res(1) : rad_a;
mesh.y = linspace(0,2*pi,N_alpha);
% mesh.y = mesh.y + mesh_res(2)/2;
mesh.z = 0 : mesh_res(3) : length;
CSX = DefineRectGrid(CSX, 1e-3,mesh);

%% fake pml %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
finalKappa = 0.3/abs_length^4;
finalSigma = finalKappa*MUE0/EPS0/epsR;
CSX = AddMaterial(CSX,'pml');
CSX = SetMaterialProperty(CSX,'pml','Kappa',finalKappa,'Epsilon',epsR);
CSX = SetMaterialProperty(CSX,'pml','Sigma',finalSigma);
CSX = SetMaterialWeight(CSX,'pml','Kappa',['pow(abs(z)-' num2str(length-abs_length) ',4)']);
CSX = SetMaterialWeight(CSX,'pml','Sigma',['pow(abs(z)-' num2str(length-abs_length) ',4)']);

start = [rad_i mesh.y(1) length-abs_length];
stop = [rad_a mesh.y(end) length];
CSX = AddBox(CSX,'pml',0 ,start,stop);


CSX = AddMaterial(CSX,'fill');
CSX = SetMaterialProperty(CSX,'fill','Epsilon',epsR);
start = [mesh.x(1) mesh.y(1) 0];
stop = [mesh.x(end) mesh.y(end) length];
CSX = AddBox(CSX,'fill',0 ,start,stop);

start = [rad_i mesh.y(1) 0];
stop = [rad_a mesh.y(end) 0];

%% apply the excitation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CSX = AddExcitation(CSX,'excite',0,[1 0 0]);
weight{1} = '1/rho';
weight{2} = 0;
weight{3} = 0;
CSX = SetExcitationWeight(CSX, 'excite', weight );
CSX = AddBox(CSX,'excite',0 ,start,stop);
 
%% define dump boxes... %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CSX = AddDump(CSX,'Et_','DumpMode',0);
start = [mesh.x(1) , 0 , mesh.z(1)];
stop = [mesh.x(end) , 0 , mesh.z(end)];
CSX = AddBox(CSX,'Et_',0 , start,stop);

CSX = AddDump(CSX,'Ht_','DumpType',1,'DumpMode',0);
CSX = AddBox(CSX,'Ht_',0,start,stop);

% voltage calc (take a voltage average to be at the same spot as the
% current calculation)
CSX = AddProbe(CSX,'ut1_1',0);
start = [ rad_i 0 port_dist ];stop = [ rad_a 0 port_dist ];
CSX = AddBox(CSX,'ut1_1', 0 ,start,stop);
CSX = AddProbe(CSX,'ut1_2',0);
start = [ rad_i 0 port_dist+mesh_res(3) ];stop = [ rad_a 0 port_dist+mesh_res(3) ];
CSX = AddBox(CSX,'ut1_2', 0 ,start,stop);

% current calc
CSX = AddProbe(CSX,'it1',1);
mid = 0.5*(rad_i+rad_a);
start = [ 0 mesh.y(1) port_dist ];stop = [ mid mesh.y(end) port_dist ];
CSX = AddBox(CSX,'it1', 0 ,start,stop);

%Write openEMS compatoble xml-file
WriteOpenEMS([Sim_Path '/' Sim_CSX],FDTD,CSX);

%% cd to working dir and run openEMS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
savePath = pwd();
cd(Sim_Path); %cd to working dir
command = [openEMS_Path 'openEMS.sh ' Sim_CSX ' ' openEMS_opts];
disp(command);
system(command)
cd(savePath);

%% postproc & do the plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
UI = ReadUI({'ut1_1','ut1_2','it1'},'tmp/');
u_f = (UI.FD{1}.val + UI.FD{2}.val)/2;          %averaging voltages to fit current
i_f = UI.FD{3}.val;

delta_t = UI.TD{3}.t(1) - UI.TD{1}.t(1);        % half time-step (s) 
i_f2 = i_f .* exp(-1i*2*pi*UI.FD{1}.f*delta_t); % compensate half time-step advance of H-field

ZL = Z0/2/pi/sqrt(epsR)*log(rad_a/rad_i); %analytic line-impedance of a coax
plot(UI.FD{1}.f,ZL*ones(size(u_f)),'g');
hold on;
grid on;
Z = u_f./i_f2;
plot(UI.FD{1}.f,real(Z),'Linewidth',2);
plot(UI.FD{1}.f,imag(Z),'r','Linewidth',2);
xlim([0 2*f0]);
legend('Z_L','\Re\{Z\}','\Im\{Z\}','Location','Best');


