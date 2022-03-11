%
% Tutorials / Dipole SAR + Power budget
%
% Description at:
% http://openems.de/index.php/Tutorial:_Dipole_SAR
%
% Tested with
%  - openEMS v0.0.33
%
% (C) 2013-2015 Thorsten Liebig <thorsten.liebig@uni-due.de>

close all
clear
clc

%% switches & options...
postprocessing_only = 0;

%% prepare simulation folder
Sim_Path = 'tmp_Dipole_SAR';
Sim_CSX = 'Dipole_SAR.xml';

%% setup the simulation
physical_constants;
unit = 1e-3; % all lengths in mm

feed.R = 50; % feed resistance

%% define phantom
phantom{1}.name='skin';
phantom{1}.epsR = 50;
phantom{1}.kappa = 0.65; % S/m
phantom{1}.density = 1100; % kg/m^3
phantom{1}.radius = [80 100 100]; % ellipsoide
phantom{1}.center = [100 0 0];

phantom{2}.name='headbone';
phantom{2}.epsR = 13;
phantom{2}.kappa = 0.1; % S/m
phantom{2}.density = 2000; % kg/m^3
phantom{2}.radius = [75 95 95]; % ellipsoide
phantom{2}.center = [100 0 0];

phantom{3}.name='brain';
phantom{3}.epsR = 60;
phantom{3}.kappa = 0.7; % S/m
phantom{3}.density = 1040; % kg/m^3
phantom{3}.radius = [65 85 85]; % ellipsoide
phantom{3}.center = [100 0 0];

%% setup FDTD parameter & excitation function
f0 = 1e9; % center frequency
lambda0 = c0/f0;

f_stop = 1.5e9; % 20 dB corner frequency
lambda_min = c0/f_stop;

mesh_res_air = lambda_min/20/unit;
mesh_res_phantom = 2.5;

dipole_length = 0.46*lambda0/unit;
disp(['Lambda-half dipole length: ' num2str(dipole_length) 'mm'])

%%
FDTD = InitFDTD();
FDTD = SetGaussExcite( FDTD, 0, f_stop );
% apply PML-8 boundary conditions in all directions
BC = {'PML_8' 'PML_8' 'PML_8' 'PML_8' 'PML_8' 'PML_8'};
FDTD = SetBoundaryCond( FDTD, BC );

%% setup CSXCAD geometry & mesh
CSX = InitCSX();

%% Dipole
CSX = AddMetal( CSX, 'Dipole' ); % create a perfect electric conductor (PEC)
CSX = AddBox(CSX, 'Dipole', 1, [0 0 -dipole_length/2], [0 0 dipole_length/2]);

% mesh lines for the dipole
mesh.x = 0;
mesh.y = 0;
mesh.z = [-dipole_length/2-[-1/3 2/3]*mesh_res_phantom dipole_length/2+[-1/3 2/3]*mesh_res_phantom];

%% add the dielectrics
for n=1:numel(phantom)
  CSX = AddMaterial( CSX, phantom{n}.name );
  CSX = SetMaterialProperty( CSX, phantom{n}.name, 'Epsilon', phantom{n}.epsR, 'Kappa', phantom{n}.kappa, 'Density', phantom{n}.density);
  CSX = AddSphere( CSX, phantom{n}.name, 10+n, [0 0 0], 1,'Transform',{'Scale',phantom{n}.radius, 'Translate', phantom{n}.center} ); 

  %% mesh lines for the dielectrics
  mesh.x = [mesh.x phantom{n}.radius(1)*[-1 1]+phantom{n}.center(1) ];
  mesh.y = [mesh.y phantom{n}.radius(2)*[-1 1]+phantom{n}.center(2) ];
  mesh.z = [mesh.z phantom{n}.radius(3)*[-1 1]+phantom{n}.center(3) ];
end

%% apply the excitation & resist as a current source
[CSX port] = AddLumpedPort(CSX, 100, 1, feed.R, [-0.1 -0.1 -mesh_res_phantom/2], [0.1 0.1 +mesh_res_phantom/2], [0 0 1], true);

% mesh lines for the port
mesh.z = [mesh.z -mesh_res_phantom/2 +mesh_res_phantom/2];

%% smooth the mesh over the dipole and phantom
mesh = SmoothMesh(mesh, mesh_res_phantom);

%% add lines for the air-box
mesh.x = [mesh.x -200 250+100];
mesh.y = [mesh.y -250 250];
mesh.z = [mesh.z -250 250];

% smooth the final mesh (incl. air box)
mesh = SmoothMesh(mesh, mesh_res_air, 1.2);

%% dump SAR
start = [-10 -100 -100];
stop = [180  100  100];
CSX = AddDump( CSX, 'SAR', 'DumpType', 20, 'Frequency', f0,'FileType',1,'DumpMode',2);
CSX = AddBox( CSX, 'SAR', 0, start, stop);

%% nf2ff calc
start = [mesh.x(1)   mesh.y(1)   mesh.z(1)];
stop  = [mesh.x(end) mesh.y(end) mesh.z(end)];
[CSX nf2ff] = CreateNF2FFBox(CSX, 'nf2ff', start, stop, 'OptResolution', lambda_min/15/unit);

%%
% add 10 equidistant cells (air)
% around the structure to keep the pml away from the nf2ff box
mesh = AddPML( mesh, 10 );

% Define the mesh
CSX = DefineRectGrid(CSX, unit, mesh);

%%
if (postprocessing_only==0)
    [status, message, messageid] = rmdir( Sim_Path, 's' ); % clear previous directory
    [status, message, messageid] = mkdir( Sim_Path ); % create empty simulation folder

    % write openEMS compatible xml-file
    WriteOpenEMS( [Sim_Path '/' Sim_CSX], FDTD, CSX );

    % show the structure
    CSXGeomPlot( [Sim_Path '/' Sim_CSX] );

    % run openEMS
    RunOpenEMS( Sim_Path, Sim_CSX );
end


%% postprocessing & make the plots
freq = linspace(500e6, 1500e6, 501 );
port = calcPort(port, Sim_Path, freq);

s11 = port.uf.ref./port.uf.inc;
Zin = port.uf.tot./port.if.tot;

Pin_f0 = interp1(freq, port.P_acc, f0);

%%
% plot feed point impedance
figure
plot( freq/1e6, real(Zin), 'k-', 'Linewidth', 2 );
hold on
grid on
plot( freq/1e6, imag(Zin), 'r--', 'Linewidth', 2 );
title( 'feed point impedance' );
xlabel( 'frequency f / MHz' );
ylabel( 'impedance Z_{in} / Ohm' );
legend( 'real', 'imag' );

% plot reflection coefficient S11
figure
plot( freq/1e9, 20*log10(abs(s11)), 'k-', 'Linewidth', 2 );
grid on
title( 'reflection coefficient' );
xlabel( 'frequency f / MHz' );
ylabel( 'S_{11} (dB)' );

%% read SAR and visualize
SAR_field = ReadHDF5Dump([Sim_Path '/SAR.h5']);

SAR = SAR_field.FD.values{1};
ptotal = ReadHDF5Attribute([Sim_Path '/SAR.h5'],'/FieldData/FD/f0','power');

%% calculate 3D pattern
phi = 0:3:360;
theta = 0:3:180;

disp( 'calculating 3D far field pattern...' );
nf2ff = CalcNF2FF(nf2ff, Sim_Path, f0, theta*pi/180, phi*pi/180, 'Outfile','3D_Pattern.h5');

%%
disp(['max SAR: ' num2str(max(SAR(:))/Pin_f0) ' W/kg normalized to 1 W accepted power']);
disp(['accepted power: ' num2str(Pin_f0) ' W (100 %)']);
disp(['radiated power: ' num2str(nf2ff.Prad) ' W ( ' num2str(round(100*(nf2ff.Prad) / Pin_f0)) ' %)']);
disp(['absorbed power: ' num2str(ptotal) ' W ( ' num2str(round(100*(ptotal) / Pin_f0)) ' %)']);
disp(['power budget:   ' num2str(100*(nf2ff.Prad + ptotal) / Pin_f0) ' %']);

%%  plot on a x/y-plane
[SAR_field SAR_mesh] = ReadHDF5Dump([Sim_Path '/SAR.h5'],'Range',{[],[],0});
figure
[X Y] = ndgrid(SAR_mesh.lines{1},SAR_mesh.lines{2});
h = pcolor(X,Y,log10(SAR_field.FD.values{1}/abs(Pin_f0)));
title( 'logarithmic SAR on an xy-plane' );
xlabel('x -->')
ylabel('y -->')
axis equal tight
set(h,'EdgeColor','none');

%%  plot on a x/z-plane
[SAR_field SAR_mesh] = ReadHDF5Dump([Sim_Path '/SAR.h5'],'Range',{[],0,[]});
figure
[X Z] = ndgrid(SAR_mesh.lines{1},SAR_mesh.lines{3});
h = pcolor(X,Z,log10(squeeze(SAR_field.FD.values{1}))/abs(Pin_f0));
title( 'logarithmic SAR on an xz-plane' );
xlabel('x -->')
ylabel('z -->')
axis equal tight
set(h,'EdgeColor','none');

%% dump SAR to vtk file
disp(['Full local/normalized SAR has been dumped to vtk file! Use Paraview to visualize']);
ConvertHDF5_VTK([Sim_Path '/SAR.h5'],[Sim_Path '/SAR'],'weight',1/abs(Pin_f0),'FieldName','SAR_local' );

