%
% Tutorials / horn antenna
%
% Describtion at:
% http://openems.de/index.php/Tutorial:_Horn_Antenna
%
% Tested with
%  - Matlab 2011a / Octave 3.6.4
%  - openEMS v0.0.31
%
% (C) 2011,2012,2013 Thorsten Liebig <thorsten.liebig@uni-due.de>

close all
clear
clc

%% setup the simulation
physical_constants;
unit = 1e-3; % all length in mm

% horn width in x-direction
horn.width  = 20;
% horn height in y-direction
horn.height = 30;
% horn length in z-direction
horn.length = 50;

horn.feed_length = 50;

horn.thickness = 2;

% horn opening angle in x, y
horn.angle = [20 20]*pi/180;

% size of the simulation box
SimBox = [200 200 200];

% frequency range of interest
f_start =  10e9;
f_stop  =  20e9;

% frequency of interest
f0 = 15e9;

%waveguide TE-mode definition
TE_mode = 'TE10';
a = horn.width;
b = horn.height;

%% setup FDTD parameter & excitation function
FDTD = InitFDTD('EndCriteria', 1e-4);
FDTD = SetGaussExcite(FDTD,0.5*(f_start+f_stop),0.5*(f_stop-f_start));
BC = {'PML_8' 'PML_8' 'PML_8' 'PML_8' 'PML_8' 'PML_8'}; % boundary conditions
FDTD = SetBoundaryCond( FDTD, BC );

%% setup CSXCAD geometry & mesh
% currently, openEMS cannot automatically generate a mesh
max_res = c0 / (f_stop) / unit / 15; % cell size: lambda/20
CSX = InitCSX();

%create fixed lines for the simulation box, substrate and port
mesh.x = [-SimBox(1)/2 -a/2 a/2 SimBox(1)/2];
mesh.x = SmoothMeshLines( mesh.x, max_res, 1.4); % create a smooth mesh between specified fixed mesh lines

mesh.y = [-SimBox(2)/2 -b/2 b/2 SimBox(2)/2];
mesh.y = SmoothMeshLines( mesh.y, max_res, 1.4 );

%create fixed lines for the simulation box and given number of lines inside the substrate
mesh.z = [-horn.feed_length 0 SimBox(3)-horn.feed_length ];
mesh.z = SmoothMeshLines( mesh.z, max_res, 1.4 );

CSX = DefineRectGrid( CSX, unit, mesh );

%% create horn
% horn feed rect waveguide
CSX = AddMetal(CSX, 'horn');
start = [-a/2-horn.thickness -b/2 mesh.z(1)];
stop  = [-a/2                 b/2 0];
CSX = AddBox(CSX,'horn',10,start,stop);
start = [a/2+horn.thickness -b/2 mesh.z(1)];
stop  = [a/2                 b/2 0];
CSX = AddBox(CSX,'horn',10,start,stop);
start = [-a/2-horn.thickness b/2+horn.thickness mesh.z(1)];
stop  = [ a/2+horn.thickness b/2                0];
CSX = AddBox(CSX,'horn',10,start,stop);
start = [-a/2-horn.thickness -b/2-horn.thickness mesh.z(1)];
stop  = [ a/2+horn.thickness -b/2                0];
CSX = AddBox(CSX,'horn',10,start,stop);

% horn opening
p(2,1) = a/2;
p(1,1) = 0;
p(2,2) = a/2 + sin(horn.angle(1))*horn.length;
p(1,2) = horn.length;
p(2,3) = -a/2 - sin(horn.angle(1))*horn.length;
p(1,3) = horn.length;
p(2,4) = -a/2;
p(1,4) = 0;
CSX = AddLinPoly( CSX, 'horn', 10, 1, -horn.thickness/2, p, horn.thickness, 'Transform', {'Rotate_X',horn.angle(2),'Translate',['0,' num2str(-b/2-horn.thickness/2) ',0']});
CSX = AddLinPoly( CSX, 'horn', 10, 1, -horn.thickness/2, p, horn.thickness, 'Transform', {'Rotate_X',-horn.angle(2),'Translate',['0,' num2str(b/2+horn.thickness/2) ',0']});

p(1,1) = b/2+horn.thickness;
p(2,1) = 0;
p(1,2) = b/2+horn.thickness + sin(horn.angle(2))*horn.length;
p(2,2) = horn.length;
p(1,3) = -b/2-horn.thickness - sin(horn.angle(2))*horn.length;
p(2,3) = horn.length;
p(1,4) = -b/2-horn.thickness;
p(2,4) = 0;
CSX = AddLinPoly( CSX, 'horn', 10, 0, -horn.thickness/2, p, horn.thickness, 'Transform', {'Rotate_Y',-horn.angle(2),'Translate',[ num2str(-a/2-horn.thickness/2) ',0,0']});
CSX = AddLinPoly( CSX, 'horn', 10, 0, -horn.thickness/2, p, horn.thickness, 'Transform', {'Rotate_Y',+horn.angle(2),'Translate',[ num2str(a/2+horn.thickness/2) ',0,0']});

% horn aperture
A = (a + 2*sin(horn.angle(1))*horn.length)*unit * (b + 2*sin(horn.angle(2))*horn.length)*unit;

%% apply the excitation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
start=[-a/2 -b/2 mesh.z(8) ];
stop =[ a/2  b/2 mesh.z(1)+horn.feed_length/2 ];
[CSX, port] = AddRectWaveGuidePort( CSX, 0, 1, start, stop, 2, a*unit, b*unit, TE_mode, 1);

%% nf2ff calc
start = [mesh.x(9) mesh.y(9) mesh.z(9)];
stop  = [mesh.x(end-8) mesh.y(end-8) mesh.z(end-8)];
[CSX nf2ff] = CreateNF2FFBox(CSX, 'nf2ff', start, stop, 'Directions', [1 1 1 1 0 1]);

%% prepare simulation folder
Sim_Path = 'tmp_Horn_Antenna';
Sim_CSX = 'horn_ant.xml';

[status, message, messageid] = rmdir( Sim_Path, 's' ); % clear previous directory
[status, message, messageid] = mkdir( Sim_Path ); % create empty simulation folder

%% write openEMS compatible xml-file
WriteOpenEMS([Sim_Path '/' Sim_CSX], FDTD, CSX);

%% show the structure
CSXGeomPlot([Sim_Path '/' Sim_CSX]);

%% run openEMS
RunOpenEMS(Sim_Path, Sim_CSX);

%% postprocessing & do the plots
freq = linspace(f_start,f_stop,201);

port = calcPort(port, Sim_Path, freq);

Zin = port.uf.tot ./ port.if.tot;
s11 = port.uf.ref ./ port.uf.inc;

plot( freq/1e9, 20*log10(abs(s11)), 'k-', 'Linewidth', 2 );
ylim([-60 0]);
grid on
title( 'reflection coefficient S_{11}' );
xlabel( 'frequency f / GHz' );
ylabel( 'reflection coefficient |S_{11}|' );

drawnow

%% NFFF contour plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate the far field at phi=0 degrees and at phi=90 degrees
thetaRange = (0:2:359) - 180;
disp( 'calculating far field at phi=[0 90] deg...' );
nf2ff = CalcNF2FF(nf2ff, Sim_Path, f0, thetaRange*pi/180, [0 90]*pi/180);

Dlog=10*log10(nf2ff.Dmax);
G_a = 4*pi*A/(c0/f0)^2;
e_a = nf2ff.Dmax/G_a;

% display some antenna parameter
disp( ['radiated power: Prad = ' num2str(nf2ff.Prad) ' Watt']);
disp( ['directivity: Dmax = ' num2str(Dlog) ' dBi'] );
disp( ['aperture efficiency: e_a = ' num2str(e_a*100) '%'] );

%%
% normalized directivity
figure
plotFFdB(nf2ff,'xaxis','theta','param',[1 2]);
drawnow
%   D_log = 20*log10(nf2ff.E_norm{1}/max(max(nf2ff.E_norm{1})));
%   D_log = D_log + 10*log10(nf2ff.Dmax);
%   plot( nf2ff.theta, D_log(:,1) ,'k-', nf2ff.theta, D_log(:,2) ,'r-' );

% polar plot
figure
polarFF(nf2ff,'xaxis','theta','param',[1 2],'logscale',[-40 20], 'xtics', 12);
drawnow
%   polar( nf2ff.theta, nf2ff.E_norm{1}(:,1) )

%% calculate 3D pattern
phiRange = sort( unique( [-180:5:-100 -100:2.5:-50 -50:1:50 50:2.5:100 100:5:180] ) );
thetaRange = sort( unique([ 0:1:50 50:2.:100 100:5:180 ]));

disp( 'calculating 3D far field...' );
nf2ff = CalcNF2FF(nf2ff, Sim_Path, f0, thetaRange*pi/180, phiRange*pi/180, 'Verbose',2,'Outfile','nf2ff_3D.h5');

figure
plotFF3D(nf2ff);

%%
E_far_normalized = nf2ff.E_norm{1}/max(nf2ff.E_norm{1}(:));
DumpFF2VTK([Sim_Path '/Horn_Pattern.vtk'],E_far_normalized,thetaRange,phiRange,'scale',1e-3);
