%
% Tutorials / radar cross section of a metal sphere
%
% Describtion at:
% http://openems.de/index.php/Tutorial:_RCS_Sphere
%
% Tested with
%  - Matlab 2011a / Octave 3.4.3
%  - openEMS v0.0.29
%
% (C) 2012 Thorsten Liebig <thorsten.liebig@uni-due.de>

close all
clear
clc

%% setup the simulation
physical_constants;
unit = 1e-3; % all length in mm

sphere.rad = 200;

% size of the simulation box
SimBox = 1000;
PW_Box = 750;

%% setup FDTD parameter & excitation function
f_start = 200e6; % start frequency
f_stop = 1000e6; % stop  frequency

FDTD = InitFDTD( 30000 );
FDTD = SetGaussExcite( FDTD, 0.5*(f_start+f_stop), 0.5*(f_stop-f_start) );
BC = [1 1 1 1 1 1]*3;  % set boundary conditions
FDTD = SetBoundaryCond( FDTD, BC );

%% setup CSXCAD geometry & mesh
% currently, openEMS cannot automatically generate a mesh
max_res = c0 / f_stop / unit / 20; % cell size: lambda/20
CSX = InitCSX();

%create fixed lines for the simulation box, substrate and port
mesh.x = SmoothMeshLines([-SimBox/2 0 SimBox/2], max_res);
mesh.y = mesh.x;
mesh.z = mesh.x;

%% create metal sphere
CSX = AddMetal( CSX, 'sphere' ); % create a perfect electric conductor (PEC)
CSX = AddSphere(CSX,'sphere',10,[0 0 0],sphere.rad);

%% plane wave excitation
k_dir = [1 0 0]; % plane wave direction --> x-direction
E_dir = [0 0 1]; % plane wave polarization --> E_z

CSX = AddPlaneWaveExcite(CSX, 'plane_wave', k_dir, E_dir);
start = [-PW_Box/2 -PW_Box/2 -PW_Box/2];
stop  = -start;
CSX = AddBox(CSX, 'plane_wave', 0, start, stop);

%% dump boxes
CSX = AddDump(CSX, 'Et');
start = [mesh.x(1)   mesh.y(1)   0];
stop  = [mesh.x(end) mesh.y(end) 0];
CSX = AddBox(CSX, 'Et', 0, start, stop);

%%nf2ff calc
start = [mesh.x(1)     mesh.y(1)     mesh.z(1)];
stop  = [mesh.x(end) mesh.y(end) mesh.z(end)];
[CSX nf2ff] = CreateNF2FFBox(CSX, 'nf2ff', start, stop);

% add 8 lines in all direction as pml spacing
mesh = AddPML(mesh,8);

CSX = DefineRectGrid( CSX, unit, mesh );

%% prepare simulation folder
Sim_Path = 'Sphere_RCS';
Sim_CSX = 'Sphere_RCS.xml';

[status, message, messageid] = rmdir( Sim_Path, 's' ); % clear previous directory
[status, message, messageid] = mkdir( Sim_Path ); % create empty simulation folder

%% write openEMS compatible xml-file
WriteOpenEMS( [Sim_Path '/' Sim_CSX], FDTD, CSX );

%% show the structure
CSXGeomPlot( [Sim_Path '/' Sim_CSX] );

%% run openEMS
RunOpenEMS( Sim_Path, Sim_CSX);

%%
disp('Use Paraview to display the elctric fields dumped by openEMS');

%%
freq = 500e6;
EF = ReadUI( 'et', Sim_Path, freq ); % time domain/freq domain voltage
Pin = 0.5*norm(E_dir)^2/Z0 .* abs(EF.FD{1}.val).^2;

%%
nf2ff = CalcNF2FF(nf2ff, Sim_Path, freq, [-180:2:180]*pi/180, [0 90]*pi/180,'Mode',1);
RCS = 4*pi./Pin(1).*nf2ff.P_rad{1}(:,1);
polar(nf2ff.theta,RCS);
xlabel('z -->');
ylabel('x -->');
hold on
grid on

%%
disp( 'calculating 3D far field pattern and dumping to vtk (use Paraview to visualize)...' );
thetaRange = (0:2:180);
phiRange = (0:2:360) - 180;
nf2ff = CalcNF2FF(nf2ff, Sim_Path, freq, thetaRange*pi/180, phiRange*pi/180,'Verbose',1,'Outfile','3D_Pattern.h5','Mode',1);

RCS = 4*pi./Pin(1).*nf2ff.P_rad{1};
DumpFF2VTK([Sim_Path '/3D_Pattern.vtk'],RCS,thetaRange,phiRange,1e-3);

disp('Use Paraview to display the 3D scattering pattern');
