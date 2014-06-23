%
% Tutorials / bi-quad antenna
%
% Tested with
%  - Octave 3.8.1
%  - openEMS v0.0.32
%
% (C) 2011-2014 Thorsten Liebig <thorsten.liebig@uni-due.de>

close all
clear
clc

%% setup the simulation
physical_constants;
unit = 1e-3; % all length in mm

quad_size = 110;
port_length = 10;
quad_mesh = 5;

Feed_R = 75;

% size of the simulation box
SimBox = [800 800 400];

% frequency range of interest
f_start =  400e6;
f_stop  =  1000e6;

% frequency of interest
f0 = 700e6;
freq = linspace(f_start,f_stop,201);

%% setup FDTD parameter & excitation function
FDTD = InitFDTD( 'endCriteria', 1e-4 );
FDTD = SetGaussExcite(FDTD,0.5*(f_start+f_stop),0.5*(f_stop-f_start));
BC = {'PML_8' 'PML_8' 'PML_8' 'PML_8' 'PML_8' 'PML_8'}; % boundary conditions
FDTD = SetBoundaryCond( FDTD, BC );

%% setup CSXCAD geometry & mesh
CSX = InitCSX();

%create fixed lines for the antenna outline and port
mesh.x = [-quad_size*sqrt(2) -quad_size/sqrt(2) 0 quad_size/sqrt(2) quad_size*sqrt(2)];
mesh.y = [-quad_size/sqrt(2)  -port_length/2 0 port_length/2 quad_size/sqrt(2)];
mesh.z = [0];

mesh = SmoothMesh(mesh, quad_mesh, 1.3);

% add air box
mesh.x = [mesh.x -SimBox(1)/2 SimBox(1)/2];
mesh.y = [mesh.y -SimBox(2)/2 SimBox(2)/2];
mesh.z = [-SimBox(3)/2 0 SimBox(3)/2];

max_res = c0 / (f_stop) / unit / 20; % cell size: lambda/20
mesh = SmoothMesh(mesh, max_res, 1.4);

CSX = DefineRectGrid( CSX, unit, mesh );

%% create bi-quad
points(1,1) = 0;
points(2,1) = port_length/2;
points(3,1) = 0;
points(1,end+1) = quad_size/sqrt(2);
points(2,end) = quad_size/sqrt(2);
points(1,end+1) = quad_size*sqrt(2);
points(2,end) = 0;
points(1,end+1) = quad_size/sqrt(2);
points(2,end) = -quad_size/sqrt(2);
points(1,end+1) = 0;
points(2,end) = -port_length/2;
points(1,end+1) = -quad_size/sqrt(2);
points(2,end) = -quad_size/sqrt(2);
points(1,end+1) = -quad_size*sqrt(2);
points(2,end) = 0;
points(1,end+1) = -quad_size/sqrt(2);
points(2,end) =  quad_size/sqrt(2);
points(1,end+1) = 0;
points(2,end) = port_length/2;

% create a thin metal wire...
CSX = AddMetal(CSX,'metal'); %create PEC with propName 'metal'
CSX = AddCurve(CSX,'metal',10, points);

%% apply the excitation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
start = [0 -port_length/2 0];
stop  = [0  port_length/2 0];
[CSX port] = AddLumpedPort(CSX,10,0,Feed_R,start,stop,[0 1 0], true);

%% nf2ff calc
start = [mesh.x(9) mesh.y(9) mesh.z(9)];
stop  = [mesh.x(end-8) mesh.y(end-8) mesh.z(end-8)];
[CSX nf2ff] = CreateNF2FFBox(CSX, 'nf2ff', start, stop);

%% prepare simulation folder
Sim_Path = 'tmp';
Sim_CSX = 'bi_quad_ant.xml';

[status, message, messageid] = rmdir( Sim_Path, 's' ); % clear previous directory
[status, message, messageid] = mkdir( Sim_Path ); % create empty simulation folder

%% write openEMS compatible xml-file
WriteOpenEMS([Sim_Path '/' Sim_CSX], FDTD, CSX);

%% show the structure
CSXGeomPlot([Sim_Path '/' Sim_CSX]);

%% run openEMS
RunOpenEMS(Sim_Path, Sim_CSX);

%% postprocessing & do the plots
port = calcPort(port, Sim_Path, freq);
s11 = port.uf.ref ./ port.uf.inc;

% plot reflection coefficient S11
figure
plot( freq/1e9, 20*log10(abs(s11)), 'k-', 'Linewidth', 2 );
ylim([-30 0]);
grid on
title( 'reflection coefficient S_{11}' );
xlabel( 'frequency f / GHz' );
ylabel( 'reflection coefficient |S_{11}|' );

%% calculate 3D far field pattern
phiRange = -180:2.5:180;
thetaRange =  0:2.5:180;

nf2ff = CalcNF2FF(nf2ff, Sim_Path, f0, thetaRange*pi/180, phiRange*pi/180);

disp( ['directivity: Dmax = ' num2str(10*log10(nf2ff.Dmax)) ' dBi'] );

% plot far-field pattern with Matlab
figure
plotFF3D(nf2ff, 'logscale', -20)

%%
disp( 'Dumping far-field pattern to vtk (use Paraview to visualize)...' );
DumpFF2VTK('Bi_Quad_Pattern.vtk', nf2ff.E_norm{1} / max(nf2ff.E_norm{1}(:)) * nf2ff.Dmax, thetaRange, phiRange, 'scale', 0.05);
