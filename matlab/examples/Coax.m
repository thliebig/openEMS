close all
clear
clc

%% setup the simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
abs_length = 250;
length = 1000;
coax_rad_i  = 100;
coax_rad_ai = 230;
coax_rad_aa = 240;
mesh_res = [5 5 5];

EPS0 = 8.85418781762e-12;
MUE0 = 1.256637062e-6;
C0 = 1/sqrt(EPS0*MUE0);
Z0 = sqrt(MUE0/EPS0);

f0 = 0.5e9;
epsR = 1;

%% define openEMS options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
openEMS_opts = '';
% openEMS_opts = [openEMS_opts ' --debug-material'];
% openEMS_opts = [openEMS_opts ' --debug-operator'];

% openEMS_opts = [openEMS_opts ' --disable-dumps --engine=fastest'];
openEMS_opts = [openEMS_opts ' --engine=sse-compressed'];

Sim_Path = 'tmp';
Sim_CSX = 'coax.xml';

if (exist(Sim_Path,'dir'))
    rmdir(Sim_Path,'s');
end
mkdir(Sim_Path);

%% setup FDTD parameter & excitation function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FDTD = InitFDTD(5e5,1e-5);
FDTD = SetGaussExcite(FDTD,f0,f0);
BC = [0 0 0 0 0 0]; %electric walls only
FDTD = SetBoundaryCond(FDTD,BC);

%% setup CSXCAD geometry & mesh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CSX = InitCSX();
mesh.x = -2.5*mesh_res(1)-coax_rad_aa : mesh_res(1) : coax_rad_aa+2.5*mesh_res(1);
mesh.y = mesh.x;
mesh.z = 0 : mesh_res(3) : length;
CSX = DefineRectGrid(CSX, 1e-3,mesh);

%% fake pml %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
finalKappa = 0.3/abs_length^4;
finalSigma = finalKappa*MUE0/EPS0;
CSX = AddMaterial(CSX,'pml');
CSX = SetMaterialProperty(CSX,'pml','Kappa',finalKappa);
CSX = SetMaterialProperty(CSX,'pml','Sigma',finalSigma);
CSX = SetMaterialWeight(CSX,'pml','Kappa',['pow(abs(z)-' num2str(length-abs_length) ',4)']);
CSX = SetMaterialWeight(CSX,'pml','Sigma',['pow(abs(z)-' num2str(length-abs_length) ',4)']);

%%% coax
CSX = AddMaterial(CSX,'copper');
CSX = SetMaterialProperty(CSX,'copper','Kappa',56e6);
start = [0, 0 , 0];stop = [0, 0 , length];
CSX = AddCylinder(CSX,'copper',0 ,start,stop,coax_rad_i);
CSX = AddCylindricalShell(CSX,'copper',0 ,start,stop,0.5*(coax_rad_aa+coax_rad_ai),(coax_rad_aa-coax_rad_ai));
start(3) = length-abs_length;
CSX = AddCylindricalShell(CSX,'pml',0 ,start,stop,0.5*(coax_rad_i+coax_rad_ai),(coax_rad_ai-coax_rad_i));

%% apply the excitation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
start(3) = 0; stop(3)=mesh_res(1)/2;
CSX = AddExcitation(CSX,'excite',0,[1 1 0]);
weight{1} = '(x)/(x*x+y*y)';
weight{2} = 'y/pow(rho,2)';
weight{3} = 0;
CSX = SetExcitationWeight(CSX, 'excite', weight );
CSX = AddCylindricalShell(CSX,'excite',0 ,start,stop,0.5*(coax_rad_i+coax_rad_ai),(coax_rad_ai-coax_rad_i));
 
%% define dump boxes... %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CSX = AddDump(CSX,'Et_','DumpMode',2);
start = [mesh.x(1) , 0 , mesh.z(1)];
stop = [mesh.x(end) , 0 , mesh.z(end)];
CSX = AddBox(CSX,'Et_',0 , start,stop);

CSX = AddDump(CSX,'Ht_','DumpType',1,'DumpMode',2);
CSX = AddBox(CSX,'Ht_',0,start,stop);

%voltage calc
CSX = AddProbe(CSX,'ut1_1',0);
start = [ coax_rad_i 0 length/2 ];stop = [ coax_rad_ai 0 length/2 ];
CSX = AddBox(CSX,'ut1_1', 0 ,start,stop);
CSX = AddProbe(CSX,'ut1_2',0);
start = [ coax_rad_i 0 length/2+mesh_res(3) ];stop = [ coax_rad_ai 0 length/2+mesh_res(3) ];
CSX = AddBox(CSX,'ut1_2', 0 ,start,stop);

%current calc
CSX = AddProbe(CSX,'it1',1);
mid = 0.5*(coax_rad_i+coax_rad_ai);
start = [ -mid -mid length/2 ];stop = [ mid mid length/2 ];
CSX = AddBox(CSX,'it1', 0 ,start,stop);

%Write openEMS compatoble xml-file
WriteOpenEMS([Sim_Path '/' Sim_CSX],FDTD,CSX);

%% cd to working dir and run openEMS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
savePath = pwd();
cd(Sim_Path); %cd to working dir
args = [Sim_CSX ' ' openEMS_opts];
invoke_openEMS(args)
cd(savePath);

%% postproc & do the plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
UI = ReadUI({'ut1_1','ut1_2','it1'},'tmp/');
u_f = (UI.FD{1}.val + UI.FD{2}.val)/2;          %averaging voltages to fit current
i_f = UI.FD{3}.val;

delta_t = UI.TD{3}.t(1) - UI.TD{1}.t(1);        % half time-step (s) 
i_f2 = i_f .* exp(-1i*2*pi*UI.FD{1}.f*delta_t); % compensate half time-step advance of H-field

ZL = Z0/2/pi/sqrt(epsR)*log(coax_rad_ai/coax_rad_i); %analytic line-impedance of a coax
plot(UI.FD{1}.f,ZL*ones(size(u_f)),'g');
hold on;
grid on;
Z = u_f./i_f2;
plot(UI.FD{1}.f,real(Z),'Linewidth',2);
plot(UI.FD{1}.f,imag(Z),'r','Linewidth',2);
xlim([0 2*f0]);
legend('Z_L','\Re\{Z\}','\Im\{Z\}','Location','Best');
