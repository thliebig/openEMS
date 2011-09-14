%
% EXAMPLE / antennas / patch antenna simple
%
% This example demonstrates how to:
%  - calculate the reflection coefficient of a patch antenna
% 
%
% Tested with
%  - Matlab 2009b
%  - Octave 3.3.52
%  - openEMS v0.0.23
%
% (C) 2010,2011 Thorsten Liebig <thorsten.liebig@uni-due.de>

close all
clear
clc

%% setup the simulation
physical_constants;
unit = 1e-3; % all length in mm

% patch width in x-direction
patch.width  = 30; % resonant length
% patch length in y-direction
patch.length = 40;

%substrate setup
substrate.epsR   = 3.38;
substrate.kappa  = 1e-3 * 2*pi*2.45e9 * EPS0*substrate.epsR;
substrate.width  = 60;
substrate.length = 60;
substrate.thickness = 1.524;
substrate.cells = 4;

%setup feeding
feed.pos = -6; %feeding position in x-direction
feed.width = 2;  %feeding port width
feed.R = 50;     %feed resistance

% size of the simulation box
SimBox = [200 200 100];

%% setup FDTD parameter & excitation function
f0 = 0e9; % center frequency
fc = 3e9; % 20 dB corner frequency (in this case 0 Hz - 3e9 Hz)
FDTD = InitFDTD( 30000, 1e-5 );
FDTD = SetGaussExcite( FDTD, f0, fc );
BC = {'MUR' 'MUR' 'MUR' 'MUR' 'MUR' 'MUR'}; % boundary conditions
FDTD = SetBoundaryCond( FDTD, BC );

%% setup CSXCAD geometry & mesh
% currently, openEMS cannot automatically generate a mesh
max_res = c0 / (f0+fc) / unit / 20; % cell size: lambda/20
CSX = InitCSX();

%create fixed lines for the simulation box, substrate and port
mesh.x = [-SimBox(1)/2 SimBox(1)/2 -substrate.width/2 substrate.width/2 -patch.width/2 patch.width/2 feed.pos];
mesh.x = SmoothMeshLines( mesh.x, max_res, 1.4); % create a smooth mesh between specified fixed mesh lines

mesh.y = [-SimBox(2)/2 SimBox(2)/2 -substrate.length/2 substrate.length/2 -feed.width/2 feed.width/2 -patch.length/2 patch.length/2];
mesh.y = SmoothMeshLines( mesh.y, max_res, 1.4 );

%create fixed lines for the simulation box and given number of lines inside the substrate
mesh.z = [-SimBox(3)/2 linspace(0,substrate.thickness,substrate.cells) SimBox(3)/2 ];
mesh.z = SmoothMeshLines( mesh.z, max_res, 1.4 );

CSX = DefineRectGrid( CSX, unit, mesh );

%% create patch
CSX = AddMetal( CSX, 'patch' ); % create a perfect electric conductor (PEC)
start = [-patch.width/2 -patch.length/2 substrate.thickness];
stop  = [ patch.width/2  patch.length/2 substrate.thickness];
CSX = AddBox(CSX,'patch',10,start,stop); % add a box-primitive to the metal property 'patch'

%% create substrate
CSX = AddMaterial( CSX, 'substrate' );
CSX = SetMaterialProperty( CSX, 'substrate', 'Epsilon', substrate.epsR, 'Kappa', substrate.kappa );
start = [-substrate.width/2 -substrate.length/2 0];
stop  = [ substrate.width/2  substrate.length/2 substrate.thickness];
CSX = AddBox( CSX, 'substrate', 0, start, stop );

%% create ground (same size as substrate)
CSX = AddMetal( CSX, 'gnd' ); % create a perfect electric conductor (PEC)
start(3)=0;
stop(3) =0;
CSX = AddBox(CSX,'gnd',10,start,stop);

%% apply the excitation & resist as a current source
start = [feed.pos-.1 -feed.width/2 0];
stop  = [feed.pos+.1 +feed.width/2 substrate.thickness];
[CSX] = AddLumpedPort(CSX, 5 ,1 ,feed.R, start, stop, [0 0 1], 'excite');

%%nf2ff calc
SimBox = SimBox - max_res * 4;  %reduced SimBox size for nf2ff box
[CSX nf2ff] = CreateNF2FFBox(CSX, 'nf2ff', -SimBox/2, SimBox/2);

%% prepare simulation folder
Sim_Path = 'tmp';
Sim_CSX = 'patch_ant.xml';

[status, message, messageid] = rmdir( Sim_Path, 's' ); % clear previous directory
[status, message, messageid] = mkdir( Sim_Path ); % create empty simulation folder

%% write openEMS compatible xml-file
WriteOpenEMS( [Sim_Path '/' Sim_CSX], FDTD, CSX );

%% show the structure
CSXGeomPlot( [Sim_Path '/' Sim_CSX] );

%% run openEMS
RunOpenEMS( Sim_Path, Sim_CSX );

%% postprocessing & do the plots
freq = linspace( max([1e9,f0-fc]), f0+fc, 501 );
U = ReadUI( {'port_ut1','et'}, 'tmp/', freq ); % time domain/freq domain voltage
I = ReadUI( 'port_it1', 'tmp/', freq ); % time domain/freq domain current (half time step is corrected)

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
if_ref = if_inc - I.FD{1}.val;
s11 = uf_ref ./ uf_inc;
plot( freq/1e6, 20*log10(abs(s11)), 'k-', 'Linewidth', 2 );
grid on
title( 'reflection coefficient S_{11}' );
xlabel( 'frequency f / MHz' );
ylabel( 'reflection coefficient |S_{11}|' );

P_in = 0.5*U.FD{1}.val .* conj( I.FD{1}.val ); % antenna feed power

%% NFFF contour plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find resonance frequncy from s11
f_res_ind = find(s11==min(s11));
f_res = freq(f_res_ind);

% calculate the far field at phi=0 degrees and at phi=90 degrees
thetaRange = (0:2:359) - 180;
r = 1; % evaluate fields at radius r
disp( 'calculating far field at phi=[0 90] deg...' );
[E_far_theta,E_far_phi,Prad,Dmax] = AnalyzeNF2FF( Sim_Path, nf2ff, f_res, thetaRange, [0 90], r );

Dlog=10*log10(Dmax);

% display power and directivity
disp( ['radiated power: Prad = ' num2str(Prad) ' Watt']);
disp( ['directivity: Dmax = ' num2str(Dlog) ' dBi'] );
disp( ['efficiency: nu_rad = ' num2str(100*Prad./real(P_in(f_res_ind))) ' %']);

% calculate the e-field magnitude for phi = 0 deg
E_phi0_far = zeros(1,numel(thetaRange));
for n=1:numel(thetaRange)
    E_phi0_far(n) = norm( [E_far_theta(n,1) E_far_phi(n,1)] );
end

E_phi0_far_log = 20*log10(abs(E_phi0_far)/max(abs(E_phi0_far)));
E_phi0_far_log = E_phi0_far_log + Dlog;

% display polar plot
figure
plot( thetaRange, E_phi0_far_log ,'k-' );
xlabel( 'theta (deg)' );
ylabel( 'directivity (dBi)');
grid on;
hold on;

% calculate the e-field magnitude for phi = 90 deg
E_phi90_far = zeros(1,numel(thetaRange));
for n=1:numel(thetaRange)
    E_phi90_far(n) = norm([E_far_theta(n,2) E_far_phi(n,2)]);
end

E_phi90_far_log = 20*log10(abs(E_phi90_far)/max(abs(E_phi90_far)));
E_phi90_far_log = E_phi90_far_log + Dlog;

% display polar plot
plot( thetaRange, E_phi90_far_log ,'r-' );
legend('phi=0','phi=90')

