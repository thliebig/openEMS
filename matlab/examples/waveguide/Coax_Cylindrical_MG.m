close all
clear
clc

%example for an cylindrical mesh, modeling a coaxial cable
% this example is using a multi-grid approach


%% setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Settings = [];
Settings.LogFile = 'openEMS.log';

physical_constants

f0 = 0.5e9;
epsR = 1;       %material filling

length = 1000;
port_dist = length/2;
rad_i  = 10;    %inner radius
rad_a = 200;    %outer radius
partial = 0.5;  %e.g. 0.5 means only one half of a coax, should be <1 or change boundary cond.
max_mesh = 10 / sqrt(epsR);
max_alpha = max_mesh;
N_alpha = ceil(rad_a * 2*pi * partial / max_alpha);

%make it even...
N_alpha = N_alpha + mod(N_alpha,2);
%make sure it is multiple of 4, needed for 2 multi-grid steps
N_alpha = ceil((N_alpha)/4) *4 + 1;

openEMS_opts = '';
% openEMS_opts = [openEMS_opts ' --disable-dumps'];
% openEMS_opts = [openEMS_opts ' --debug-material'];
% openEMS_opts = [openEMS_opts ' --numThreads=1'];

def_refSimu = 0; % do a reference simulation without the multi-grid

%% setup done %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if (def_refSimu>0)
    Sim_Path = 'tmp_ref';
else
    Sim_Path = 'tmp';
end
Sim_CSX = 'coax.xml';

if (exist(Sim_Path,'dir'))
    rmdir(Sim_Path,'s');
end
mkdir(Sim_Path);

%setup FDTD parameter
if (def_refSimu>0)
    FDTD = InitCylindricalFDTD(1e5,1e-5,'OverSampling',5 );
else
    FDTD = InitCylindricalFDTD(1e5,1e-5,'OverSampling',5 ,'MultiGrid','60,120');
end
FDTD = SetGaussExcite(FDTD,f0,f0);
BC = [0 0 1 1 2 2];
FDTD = SetBoundaryCond(FDTD,BC);

mesh_res = [max_mesh 2*pi*partial/N_alpha max_mesh];

%setup CSXCAD geometry
CSX = InitCSX();
mesh.x = SmoothMeshLines([rad_i rad_a],mesh_res(1));
mesh.y = linspace(-pi*partial,pi*partial,N_alpha);
mesh.z = SmoothMeshLines([0 port_dist length],mesh_res(3));
CSX = DefineRectGrid(CSX, 1e-3,mesh);

start = [rad_i mesh.y(1)   mesh.z(3)];
stop  = [rad_a mesh.y(end) mesh.z(3)];

CSX = AddExcitation(CSX,'excite',0,[1 0 0]);
weight{1} = '1/rho';
weight{2} = 0;
weight{3} = 0;
CSX = SetExcitationWeight(CSX, 'excite', weight );
CSX = AddBox(CSX,'excite',0 ,start,stop);
 

start = [mesh.x(1)   mesh.y(1)   mesh.z(1)];
stop  = [mesh.x(end) mesh.y(end) mesh.z(end)];
CSX = AddMaterial(CSX,'material');
CSX = SetMaterialProperty(CSX,'material','Epsilon',epsR);
CSX = AddBox(CSX,'material',0 ,start,stop);

%dump
CSX = AddDump(CSX,'Et_rz_','DumpMode',0);
start = [mesh.x(1)   0 mesh.z(1)];
stop  = [mesh.x(end) 0 mesh.z(end)];
CSX = AddBox(CSX,'Et_rz_',0 , start,stop);

CSX = AddDump(CSX,'Ht_rz_','DumpType',1,'DumpMode',0);
CSX = AddBox(CSX,'Ht_rz_',0 , start,stop);

CSX = AddDump(CSX,'Et_','DumpType',0,'DumpMode',0);
start = [mesh.x(1)   mesh.y(1)   length/2];
stop  = [mesh.x(end) mesh.y(end) length/2];
CSX = AddBox(CSX,'Et_',0,start,stop);

CSX = AddDump(CSX,'Ht_','DumpType',1,'DumpMode',0);
start = [mesh.x(1)   mesh.y(1)   length/2];
stop  = [mesh.x(end) mesh.y(end) length/2];
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
mid = 75;
start = [ 0 mesh.y(1) port_dist+mesh_res(3)/2 ];stop = [ mid mesh.y(end) port_dist+mesh_res(3)/2 ];
CSX = AddBox(CSX,'it1', 0 ,start,stop);

%% Write openEMS compatoble xml-file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WriteOpenEMS([Sim_Path '/' Sim_CSX],FDTD,CSX);

RunOpenEMS(Sim_Path, Sim_CSX, openEMS_opts, Settings)

%%
close all
freq = linspace(0,2*f0,201);
UI = ReadUI({'ut1_1','ut1_2','it1'},Sim_Path,freq);
u_f = (UI.FD{1}.val + UI.FD{2}.val)/2;          %averaging voltages to fit current
i_f = UI.FD{3}.val / partial;

% plot(UI.TD{1}.t,UI.TD{1}.val);
% grid on;
% 
% figure
% plot(UI.TD{3}.t,UI.TD{3}.val);
% grid on;

%plot Z_L compare
figure
ZL = Z0/2/pi/sqrt(epsR)*log(rad_a/rad_i); %analytic line-impedance of a coax
plot(UI.FD{1}.f,ZL*ones(size(u_f)),'g','Linewidth',3);
hold on;
grid on;
Z = u_f./i_f;
plot(UI.FD{1}.f,real(Z),'k--','Linewidth',2);
plot(UI.FD{1}.f,imag(Z),'r-','Linewidth',2);
xlim([0 2*f0]);
legend('Z_L - analytic','\Re\{Z\} - FDTD','\Im\{Z\} - FDTD','Location','Best');


