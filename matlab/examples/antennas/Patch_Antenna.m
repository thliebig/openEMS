%
% EXAMPLE / antennas / patch antenna
%
% This example demonstrates how to:
%  - calculate the reflection coefficient of a patch antenna
% 
%
% Tested with
%  - Matlab 2009b
%  - Octave 3.3.52
%  - openEMS v0.0.14
%
% (C) 2010 Thorsten Liebig <thorsten.liebig@uni-due.de>

close all
clear
clc

%% setup the simulation
physical_constants;
unit = 1e-3; % all length in mm

% width in x-direction
% length in y-direction
% main radiation in z-direction
patch.width  = 32.86; % resonant length
patch.length = 41.37;

substrate.epsR   = 3.38;
substrate.width  = 120;
substrate.length = 120;
substrate.thickness = 1.524;
substrate.cells = 5;

feed.pos = -4.5;
feed.width = 0.5;
feed.R = 50; % feed resistance

% size of the simulation box
SimBox = [120 120 32];

%% prepare simulation folder
Sim_Path = 'tmp';
Sim_CSX = 'patch_ant.xml';
[status, message, messageid] = rmdir( Sim_Path, 's' ); % clear previous directory
[status, message, messageid] = mkdir( Sim_Path ); % create empty simulation folder

%% setup FDTD parameter & excitation function
max_timesteps = 30000;
min_decrement = 1e-5; % equivalent to -50 dB
f0 = 2e9; % center frequency
fc = 1e9; % 10 dB corner frequency (in this case 1e9 Hz - 3e9 Hz)
FDTD = InitFDTD( max_timesteps, min_decrement );
FDTD = SetGaussExcite( FDTD, f0, fc );
BC = {'MUR' 'MUR' 'MUR' 'MUR' 'PEC' 'MUR'}; % boundary conditions
% BC = {'PML_8' 'PML_8' 'PML_8' 'PML_8' 'PEC' 'PML_8'}; % use pml instead of mur
FDTD = SetBoundaryCond( FDTD, BC );

%% setup CSXCAD geometry & mesh
% currently, openEMS cannot automatically generate a mesh
max_res = c0 / (f0+fc) / unit / 20; % cell size: lambda/20
CSX = InitCSX();
mesh.x = [-SimBox(1)/2 SimBox(1)/2 -substrate.width/2 substrate.width/2 feed.pos-feed.width/2 feed.pos+feed.width/2];
% add patch mesh with 2/3 - 1/3 rule
mesh.x = [mesh.x -patch.width/2-max_res*0.66 -patch.width/2+max_res*0.33 patch.width/2+max_res*0.66 patch.width/2-max_res*0.33];
mesh.x = SmoothMeshLines( mesh.x, max_res ); % create a smooth mesh between specified mesh lines
mesh.y = [-SimBox(2)/2 SimBox(2)/2 -substrate.length/2 substrate.length/2 -feed.width/2 feed.width/2];
% add patch mesh with 2/3 - 1/3 rule
mesh.y = [mesh.y -patch.length/2-max_res*0.66 -patch.length/2+max_res*0.33 patch.length/2+max_res*0.66 patch.length/2-max_res*0.33];
mesh.y = SmoothMeshLines( mesh.y, max_res );
mesh.z = [linspace(0,substrate.thickness,substrate.cells) SimBox(3) SimBox(3)];
mesh.z = SmoothMeshLines( mesh.z, max_res );
mesh = AddPML( mesh, [8 8 8 8 0 8] ); % add equidistant cells (air around the structure)
CSX = DefineRectGrid( CSX, unit, mesh );

%% create patch
CSX = AddMetal( CSX, 'patch' ); % create a perfect electric conductor (PEC)
start = [-patch.width/2 -patch.length/2 substrate.thickness];
stop  = [ patch.width/2  patch.length/2 substrate.thickness];
CSX = AddBox(CSX,'patch',10,start,stop);

%% create substrate
CSX = AddMaterial( CSX, 'substrate' );
CSX = SetMaterialProperty( CSX, 'substrate', 'Epsilon', substrate.epsR );
start = [-substrate.width/2 -substrate.length/2 0];
stop  = [ substrate.width/2  substrate.length/2 substrate.thickness];
CSX = AddBox( CSX, 'substrate', 0, start, stop );

%% apply the excitation & resist as a current source
% this creates a "port"
CSX = AddMaterial( CSX, 'resist' );
kappa = substrate.thickness/feed.R/feed.width^2/unit;
CSX = SetMaterialProperty( CSX, 'resist', 'Kappa', kappa );
start = [feed.pos-feed.width/2 -feed.width/2 0];
stop  = [feed.pos+feed.width/2  feed.width/2 substrate.thickness];
CSX = AddBox( CSX, 'resist', 15, start, stop );

CSX = AddExcitation( CSX, 'excite', 0, [0 0 1] ); % excitation in z-direction
CSX = AddBox( CSX, 'excite', 0, start, stop );

%% define voltage calc boxes
CSX = AddProbe( CSX, 'ut1', 0 );
start = [feed.pos 0 0];
stop  = [feed.pos 0 substrate.thickness];
CSX = AddBox( CSX, 'ut1', 0 , stop, start );

%% define current calc boxes
CSX = AddProbe( CSX, 'it1', 1 );
start = [feed.pos-feed.width -feed.width substrate.thickness/2];
stop  = [feed.pos+feed.width  feed.width substrate.thickness/2];
CSX = AddBox( CSX, 'it1', 0, start, stop );

%% dump magnetic field over the patch antenna
CSX = AddDump( CSX, 'Ht_', 'DumpType', 1, 'DumpMode', 2 ); % cell interpolated
start = [-patch.width -patch.length substrate.thickness+1];
stop  = [ patch.width  patch.length substrate.thickness+1];
CSX = AddBox( CSX, 'Ht_', 0, start, stop );

%% write openEMS compatible xml-file
WriteOpenEMS( [Sim_Path '/' Sim_CSX], FDTD, CSX );

%% show the structure
CSXGeomPlot( [Sim_Path '/' Sim_CSX] );

%% run openEMS
openEMS_opts = '';
openEMS_opts = [openEMS_opts ' --engine=fastest'];
RunOpenEMS( Sim_Path, Sim_CSX, openEMS_opts );

%% postprocessing & do the plots
freq = linspace( f0-fc, f0+fc, 501 );
U = ReadUI( {'ut1','et'}, 'tmp/', freq ); % time domain/freq domain voltage
I = ReadUI( 'it1', 'tmp/', freq ); % time domain/freq domain current (half time step is corrected)

% plot time domain voltage
figure
[ax,h1,h2] = plotyy( U.TD{1}.t/1e-9, U.TD{1}.val, U.TD{2}.t/1e-9, U.TD{2}.val );
set( h1, 'Linewidth', 2 );
set( h1, 'Color', [1 0 0] );
set( h2, 'Linewidth', 2 );
set( h2, 'Color', [0 0 0] );
grid on
title( 'time domain voltage' );
xlabel( 'time t / ns' );
ylabel( ax(1), 'voltage ut1 / V' );
ylabel( ax(2), 'voltage et / V' );
% now make the y-axis symmetric to y=0 (align zeros of y1 and y2)
y1 = ylim(ax(1));
y2 = ylim(ax(2));
ylim( ax(1), [-max(abs(y1)) max(abs(y1))] );
ylim( ax(2), [-max(abs(y2)) max(abs(y2))] );

% plot feed point impedance
figure
Zin = U.FD{1}.val ./ I.FD{1}.val;
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
uf_inc = 0.5*(U.FD{1}.val + I.FD{1}.val * 50);
if_inc = 0.5*(I.FD{1}.val - U.FD{1}.val / 50);
uf_ref = U.FD{1}.val - uf_inc;
if_ref = I.FD{1}.val - if_inc;
s11 = uf_ref ./ uf_inc;
plot( freq/1e6, 20*log10(abs(s11)), 'k-', 'Linewidth', 2 );
grid on
title( 'reflection coefficient S_{11}' );
xlabel( 'frequency f / MHz' );
ylabel( 'reflection coefficient |S_{11}|' );

%% visualize magnetic fields
% you will find vtk dump files in the simulation folder (tmp/)
% use paraview to visulaize them
