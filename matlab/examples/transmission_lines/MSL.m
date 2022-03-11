%
% EXAMPLE / microstrip / MSL
%
% Microstrip line on air "substrate" in z-direction.
%
% This example demonstrates:
%  - simple microstrip geometry
%  - characteristic impedance
%  - material grading function
%  - geometric priority concept
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

% geometry
abs_length = 100; % absorber length
length = 600;
width = 400;
height = 200;
MSL_width = 50;
MSL_height = 10;

%% prepare simulation folder
Sim_Path = 'tmp';
Sim_CSX = 'msl.xml';
[status, message, messageid] = rmdir( Sim_Path, 's' ); % clear previous directory
[status, message, messageid] = mkdir( Sim_Path ); % create empty simulation folder

%% setup FDTD parameter & excitation function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_timesteps = 2000;
min_decrement = 1e-5; % equivalent to -50 dB
f0 = 2e9; % center frequency
fc = 1e9; % 10 dB corner frequency (in this case 1e9 Hz - 3e9 Hz)
FDTD = InitFDTD( max_timesteps, min_decrement );
FDTD = SetGaussExcite( FDTD, f0, fc );
BC = {'PMC' 'PMC' 'PEC' 'PMC' 'PEC' 'PEC'};
FDTD = SetBoundaryCond( FDTD, BC );

%% setup CSXCAD geometry & mesh
% very simple mesh
CSX = InitCSX();
resolution = c0/(f0+fc)/unit /15; % resolution of lambda/15
mesh.x = SmoothMeshLines( [-width/2, width/2, -MSL_width/2, MSL_width/2], resolution ); % create smooth lines from fixed lines
mesh.y = SmoothMeshLines( [linspace(0,MSL_height,5) MSL_height+1 height], resolution );
mesh.z = SmoothMeshLines( [0 length], resolution );
CSX = DefineRectGrid( CSX, unit, mesh );

%% create MSL
% attention! the skin effect is not simulated, because the MSL is
% discretized with only one cell!
CSX = AddMaterial( CSX, 'copper' );
CSX = SetMaterialProperty( CSX, 'copper', 'Kappa', 56e6 );
start = [-MSL_width/2, MSL_height,   0];
stop  = [ MSL_width/2, MSL_height+1, length];
priority = 100; % the geometric priority is set to 100
CSX = AddBox( CSX, 'copper', priority, start, stop );

%% add excitation below the strip
start = [-MSL_width/2, 0         , mesh.z(1)];
stop  = [ MSL_width/2, MSL_height, mesh.z(1)];
CSX = AddExcitation( CSX, 'excite', 0, [0 -1 0] );
CSX = AddBox( CSX, 'excite', 0, start, stop );

%% fake pml
% this "pml" is a normal material with graded losses
% electric and magnetic losses are related to give low reflection
% for normally incident TEM waves
finalKappa = 1/abs_length^2;
finalSigma = finalKappa*MUE0/EPS0;
CSX = AddMaterial( CSX, 'fakepml' );
CSX = SetMaterialProperty( CSX, 'fakepml', 'Kappa', finalKappa );
CSX = SetMaterialProperty( CSX, 'fakepml', 'Sigma', finalSigma );
CSX = SetMaterialWeight( CSX, 'fakepml', 'Kappa', ['pow(z-' num2str(length-abs_length) ',2)'] );
CSX = SetMaterialWeight( CSX, 'fakepml', 'Sigma', ['pow(z-' num2str(length-abs_length) ',2)'] );
start = [mesh.x(1)  mesh.y(1)  length-abs_length];
stop  = [mesh.x(end) mesh.y(end) length];
% the geometric priority is set to 0, which is lower than the priority
% of the MSL, thus the MSL (copper) has precedence
priority = 0;
CSX = AddBox( CSX, 'fakepml', priority, start, stop );

%% define dump boxes
start = [mesh.x(1),  MSL_height/2, mesh.z(1)];
stop  = [mesh.x(end), MSL_height/2, mesh.z(end)];
CSX = AddDump( CSX, 'Et_', 'DumpMode', 2 ); % cell interpolated
CSX = AddBox( CSX, 'Et_', 0, start, stop );
CSX = AddDump( CSX, 'Ht_', 'DumpType', 1, 'DumpMode', 2 ); % cell interpolated
CSX = AddBox( CSX, 'Ht_', 0, start, stop );

%% define voltage calc box
% voltage calc boxes will automatically snap to the next mesh-line
CSX = AddProbe( CSX, 'ut1', 0 );
zidx  = interp1( mesh.z, 1:numel(mesh.z), length/2, 'nearest' );
start = [0 MSL_height mesh.z(zidx)];
stop  = [0 0          mesh.z(zidx)];
CSX = AddBox( CSX, 'ut1', 0, start, stop );
% add a second voltage probe to compensate space offset between voltage and
% current
CSX = AddProbe( CSX, 'ut2', 0 );
start = [0 MSL_height mesh.z(zidx+1)];
stop  = [0 0          mesh.z(zidx+1)];
CSX = AddBox( CSX, 'ut2', 0, start, stop );

%% define current calc box
% current calc boxes will automatically snap to the next dual mesh-line
CSX = AddProbe( CSX, 'it1', 1 );
xidx1  = interp1( mesh.x, 1:numel(mesh.x), -MSL_width/2, 'nearest' );
xidx2  = interp1( mesh.x, 1:numel(mesh.x),  MSL_width/2, 'nearest' );
xdelta = diff(mesh.x);
yidx1  = interp1( mesh.y, 1:numel(mesh.y), MSL_height, 'nearest' );
yidx2  = interp1( mesh.y, 1:numel(mesh.y), MSL_height+1, 'nearest' );
ydelta = diff(mesh.y);
zdelta = diff(mesh.z);
start = [mesh.x(xidx1)-xdelta(xidx1-1)/2, mesh.y(yidx1)-ydelta(yidx1-1)/2, mesh.z(zidx)+zdelta(zidx)/2];
stop  = [mesh.x(xidx2)+xdelta(xidx2)/2,   mesh.y(yidx2)+ydelta(yidx2)/2,   mesh.z(zidx)+zdelta(zidx)/2];
CSX = AddBox( CSX, 'it1', 0, start, stop );

%% write openEMS compatible xml-file
WriteOpenEMS( [Sim_Path '/' Sim_CSX], FDTD, CSX );

%% show the structure
CSXGeomPlot( [Sim_Path '/' Sim_CSX] );

%% run openEMS
openEMS_opts = '';
openEMS_opts = [openEMS_opts ' --engine=fastest'];
% openEMS_opts = [openEMS_opts ' --debug-material'];
% openEMS_opts = [openEMS_opts ' --debug-boxes'];
RunOpenEMS( Sim_Path, Sim_CSX, openEMS_opts );

%% postprocess
freq = linspace( f0-fc, f0+fc, 501 );
U = ReadUI( {'ut1','ut2','et'}, 'tmp/', freq ); % time domain/freq domain voltage
I = ReadUI( 'it1', 'tmp/', freq ); % time domain/freq domain current (half time step offset is corrected)

% plot time domain voltage
figure
[ax,h1,h2] = plotyy( U.TD{1}.t/1e-9, U.TD{1}.val, U.TD{3}.t/1e-9, U.TD{3}.val );
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

% calculate characteristic impedance
% arithmetic mean of ut1 and ut2 -> voltage in the middle of ut1 and ut2
U = (U.FD{1}.val + U.FD{2}.val) / 2;
Z = U ./ I.FD{1}.val;

% plot characteristic impedance
figure
plot( freq/1e6, real(Z), 'k-', 'Linewidth', 2 );
hold on
grid on
plot( freq/1e6, imag(Z), 'r--', 'Linewidth', 2 );
title( 'characteristic impedance of MSL' );
xlabel( 'frequency f / MHz' );
ylabel( 'characteristic impedance Z / Ohm' );
legend( 'real', 'imag' );

%% visualize electric and magnetic fields
% you will find vtk dump files in the simulation folder (tmp/)
% use paraview to visualize them
