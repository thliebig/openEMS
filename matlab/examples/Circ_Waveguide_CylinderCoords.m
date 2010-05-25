close all
clear
clc

%% setup the simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
abs_length = 500;
length = 10000;
unit = 1e-3;
rad  = 300;
mesh_max = 15;
N_alpha = ceil(rad * pi / mesh_max) * 2;
mesh_res = [mesh_max 2*pi/N_alpha mesh_max];

do_Half_Waveguide = 1;

EPS0 = 8.85418781762e-12;
MUE0 = 1.256637062e-6;
C0 = 1/sqrt(EPS0*MUE0);

f0 = 400e6;

p11 = 1.841;
kc = p11 / rad /unit;
k = 2*pi*f0/C0;
fc = C0*kc/2/pi;
beta = sqrt(k^2 - kc^2);

kc = kc*unit;
func_Er = [ num2str(-1/kc^2,15) '/rho*cos(a)*j1('  num2str(kc,15) '*rho)'];
func_Ea = [ num2str(1/kc,15) '*sin(a)*0.5*(j0('  num2str(kc,15) '*rho)-jn(2,'  num2str(kc,15) '*rho))'];

%% define openEMS options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
openEMS_opts = '';
% openEMS_opts = [openEMS_opts ' --disable-dumps'];
% openEMS_opts = [openEMS_opts ' --debug-material'];
% openEMS_opts = [openEMS_opts ' --debug-operator'];

if (do_Half_Waveguide)
    Sim_Path = 'tmp_half_CWG_CC';
else
    Sim_Path = 'tmp_full_CWG_CC';
end
Sim_CSX = 'Circ_WG_CC.xml';

if (exist(Sim_Path,'dir'))
    rmdir(Sim_Path,'s');
end
mkdir(Sim_Path);

%% setup FDTD parameter & excitation function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FDTD = InitCylindricalFDTD(1e5,1e-5,'OverSampling',10);
% T = 1/f0;
% FDTD = SetCustomExcite(FDTD,f0,[ '(1-exp(-1*(t/' num2str(T) ')^2) ) * sin(2*pi*' num2str(f0) '*t)' ]);
FDTD = SetSinusExcite(FDTD,f0);
BC = [0 0 0 0 0 0];
FDTD = SetBoundaryCond(FDTD,BC);

%% setup CSXCAD geometry & mesh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CSX = InitCSX();
mesh.x = 0:mesh_res(1):rad;
if (do_Half_Waveguide)
    mesh.y = linspace(-pi/2,pi/2,N_alpha/2);
else
    mesh.y = linspace(-pi,pi,N_alpha)+pi/2;
end
y_delta = mesh.y(2) - mesh.y(1);
mesh.z = 0 : mesh_res(3) : length;
CSX = DefineRectGrid(CSX, 1e-3,mesh);

%% fake pml %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
start = [0 mesh.y(1)-y_delta length-abs_length];
stop = [rad*1.2 mesh.y(end)+y_delta length];
finalKappa = 0.3/abs_length^4;
finalSigma = finalKappa*MUE0/EPS0;
CSX = AddMaterial(CSX,'pml');
CSX = SetMaterialProperty(CSX,'pml','Kappa',finalKappa);
CSX = SetMaterialProperty(CSX,'pml','Sigma',finalSigma);
CSX = SetMaterialWeight(CSX,'pml','Kappa',['pow(abs(z)-' num2str(length-abs_length) ',4)']);
CSX = SetMaterialWeight(CSX,'pml','Sigma',['pow(abs(z)-' num2str(length-abs_length) ',4)']);
CSX = AddBox(CSX,'pml',0 ,start,stop);

%% apply the excitation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CSX = AddExcitation(CSX,'excite',0,[1 1 0]);
weight{1} = func_Er;
weight{2} = func_Ea;
weight{3} = 0;
CSX = SetExcitationWeight(CSX, 'excite', weight );
start(3)=-.5;
stop(3)=0.5;
CSX = AddBox(CSX,'excite', 5 ,start,stop);
 
%% define dump boxes... %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CSX = AddDump(CSX,'Et','FileType',1,'DumpMode',0,'SubSampling','2,2,5');
start = [mesh.x(1) , mesh.y(1)-y_delta , 0];
stop = [mesh.x(end) , mesh.y(end)+y_delta , length];
CSX = AddBox(CSX,'Et',0 , start,stop);

CSX = AddDump(CSX,'Ht','FileType',1,'DumpType',1,'DumpMode',0,'SubSampling','2,2,5');
CSX = AddBox(CSX,'Ht',0 , start,stop);

% % dumpt r-z-plane to vtk-file
% CSX = AddDump(CSX,'Et_rz_','FileType',0,'DumpMode',2,'SubSampling','1,1,5');
% start = [mesh.x(1) , 0 , 0];
% stop = [mesh.x(end) , 0 , length];
% CSX = AddBox(CSX,'Et_rz_',0 , start,stop);

%% define voltage calc boxes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CSX = AddProbe(CSX,'ut_exc',0);
start = [ 0 0 0 ];stop = [ rad 0 0 ];
CSX = AddBox(CSX,'ut_exc', 0 ,start,stop);

CSX = AddProbe(CSX,'ut_1',0);
start = [ 0 0 length/2 ];stop = [ rad 0 length/2 ];
CSX = AddBox(CSX,'ut_1', 0 ,start,stop);

%% Write openEMS compatoble xml-file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WriteOpenEMS([Sim_Path '/' Sim_CSX],FDTD,CSX);

%% cd to working dir and run openEMS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
savePath = pwd();
cd(Sim_Path); %cd to working dir
args = [Sim_CSX ' ' openEMS_opts];
invoke_openEMS(args);
cd(savePath);

%% do the plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
UI = ReadUI('ut_1',[Sim_Path '/']);
plot(UI.TD{1}.t,UI.TD{1}.val);
grid on;

file = [Sim_Path '/Et.h5'];
z_planes = 1;
timestep = 10;
for z =z_planes
    figure
    if exist(file,'file')
        mesh = ReadHDF5Mesh(file);
        fields = ReadHDF5FieldData(file);

        [ALPHA RHO] = meshgrid(double(mesh.lines{1}),double(mesh.lines{2}));
        X = RHO.*cos(ALPHA);
        Y = RHO.*sin(ALPHA);

        Er = double( fields.values{timestep}(:,:,z,1) );
        Ea = double( fields.values{timestep}(:,:,z,2) );
        Ez = double( fields.values{timestep}(:,:,z,3) );

        Ex = Er.*cos(ALPHA) - Ea.*sin(ALPHA);
        Ey = Er.*sin(ALPHA) + Ea.*cos(ALPHA);

        quiver(X,Y,Ex,Ey)
        axis equal
        title(['z : ' num2str(mesh.lines{2}(z)) ' ts: ' int2str(n)] );
        Ex(10,5)
        pause(1)
    end
end
