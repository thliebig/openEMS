%
% Tutorials / conical horn antenna
%
% Describtion at:
% http://openems.de/index.php/Tutorial:_Conical_Horn_Antenna
%
% Tested with
%  - Matlab 2011a / Octave 3.4.3
%  - openEMS v0.0.27
%
% (C) 2011,2012 Thorsten Liebig <thorsten.liebig@uni-due.de>

close all
clear
clc

%% setup the simulation
physical_constants;
unit = 1e-3; % all length in mm

% horn radius
horn.radius  = 20;
% horn length in z-direction
horn.length = 50;

horn.feed_length = 50;

horn.thickness = 2;

% horn opening angle
horn.angle = 20*pi/180;

% size of the simulation box
SimBox = [100 100 100]*2;

% frequency range of interest
f_start =  10e9;
f_stop  =  20e9;

% frequency of interest
f0 = 15e9;

%% mode functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% by David M. Pozar, Microwave Engineering, third edition, page 113
freq = linspace(f_start,f_stop,201);

p11 = 1.841;
kc = p11 / horn.radius /unit;
k = 2*pi*freq/C0;
fc = C0*kc/2/pi;
beta = sqrt(k.^2 - kc^2);
ZL_a = k * Z0 ./ beta;    %analytic waveguide impedance

% mode profile E- and H-field
kc = kc*unit;
func_Er = [ num2str(-1/kc^2,'%14.13f') '/rho*cos(a)*j1('  num2str(kc,'%14.13f') '*rho)'];
func_Ea = [ num2str(1/kc,'%14.13f') '*sin(a)*0.5*(j0('  num2str(kc,'%14.13f') '*rho)-jn(2,'  num2str(kc,'%14.13f') '*rho))'];
func_Ex = ['(' func_Er '*cos(a) - ' func_Ea '*sin(a) ) * (rho<' num2str(horn.radius) ')'];
func_Ey = ['(' func_Er '*sin(a) + ' func_Ea '*cos(a) ) * (rho<' num2str(horn.radius) ')'];

func_Ha = [ num2str(-1/kc^2,'%14.13f') '/rho*cos(a)*j1('  num2str(kc,'%14.13f') '*rho)'];
func_Hr = [ '-1*' num2str(1/kc,'%14.13f') '*sin(a)*0.5*(j0('  num2str(kc,'%14.13f') '*rho)-jn(2,'  num2str(kc,'%14.13f') '*rho))'];
func_Hx = ['(' func_Hr '*cos(a) - ' func_Ha '*sin(a) ) * (rho<' num2str(horn.radius) ')'];
func_Hy = ['(' func_Hr '*sin(a) + ' func_Ha '*cos(a) ) * (rho<' num2str(horn.radius) ')'];

disp([' Cutoff frequencies for this mode and wavguide is: ' num2str(fc/1e9) ' GHz']);

if (f_start<fc)
    warning('openEMS:example','f_start is smaller than the cutoff-frequency, this may result in a long simulation... ');
end

%% setup FDTD parameter & excitation function
FDTD = InitFDTD( 30000, 1e-4 );
FDTD = SetGaussExcite(FDTD,0.5*(f_start+f_stop),0.5*(f_stop-f_start));
BC = {'PML_8' 'PML_8' 'PML_8' 'PML_8' 'PML_8' 'PML_8'}; % boundary conditions
FDTD = SetBoundaryCond( FDTD, BC );

%% setup CSXCAD geometry & mesh
% currently, openEMS cannot automatically generate a mesh
max_res = c0 / (f_stop) / unit / 15; % cell size: lambda/20
CSX = InitCSX();

%create fixed lines for the simulation box, substrate and port
mesh.x = [-SimBox(1)/2 -horn.radius 0 horn.radius SimBox(1)/2];
mesh.x = SmoothMeshLines( mesh.x, max_res, 1.4); % create a smooth mesh between specified fixed mesh lines

mesh.y = mesh.x;

%create fixed lines for the simulation box and given number of lines inside the substrate
mesh.z = [-horn.feed_length 0 SimBox(3) ];
mesh.z = SmoothMeshLines( mesh.z, max_res, 1.4 );

CSX = DefineRectGrid( CSX, unit, mesh );

%% create horn
% horn + waveguide, defined by a rotational polygon
CSX = AddMetal(CSX, 'Conical_Horn');
p(1,1) = horn.radius+horn.thickness;   % x-coord point 1
p(2,1) = -horn.feed_length;     % z-coord point 1
p(1,end+1) = horn.radius+horn.thickness;   % x-coord point 1
p(2,end) = 0;     % z-coord point 1
p(1,end+1) = horn.radius+horn.thickness + sin(horn.angle)*horn.length; % x-coord point 2
p(2,end) = horn.length; % y-coord point 2
p(1,end+1) = horn.radius + sin(horn.angle)*horn.length; % x-coord point 2
p(2,end) = horn.length; % y-coord point 2
p(1,end+1) = horn.radius;  % x-coord point 1
p(2,end) = 0;     % z-coord point 1
p(1,end+1) = horn.radius;   % x-coord point 1
p(2,end) = -horn.feed_length;     % z-coord point 1
CSX = AddRotPoly(CSX,'Conical_Horn',10,0,2,p);

% horn aperture
A = pi*((horn.radius + sin(horn.angle)*horn.length)*unit)^2;

% %% apply the excitation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% xy-mode profile excitation located directly on top of pml (first 8 z-lines)
CSX = AddExcitation(CSX,'excite',0,[1 1 0]);
weight{1} = func_Ex;
weight{2} = func_Ey;
weight{3} = 0;
CSX = SetExcitationWeight(CSX,'excite',weight);
start=[0 0 mesh.z(8)-0.1 ];
stop =[0 0 mesh.z(8)+0.1 ];
CSX = AddCylinder(CSX,'excite',0 ,start,stop,horn.radius);

CSX = AddDump(CSX,'Exc_dump');
start=[-horn.radius -horn.radius mesh.z(8)-0.1 ];
stop =[+horn.radius +horn.radius mesh.z(8)+0.1 ];
CSX = AddBox(CSX,'Exc_dump',0,start,stop);

%% voltage and current definitions using the mode matching probes %%%%%%%%%
%port 1
start = [-horn.radius -horn.radius mesh.z(1)+horn.feed_length/2];
stop  = [ horn.radius  horn.radius mesh.z(1)+horn.feed_length/2];
CSX = AddProbe(CSX, 'ut1', 10, 1, [], 'ModeFunction',{func_Ex,func_Ey,0});
CSX = AddBox(CSX,  'ut1',  0 ,start,stop);
CSX = AddProbe(CSX,'it1', 11, 1, [], 'ModeFunction',{func_Hx,func_Hy,0});
CSX = AddBox(CSX,'it1', 0 ,start,stop);

%% nf2ff calc
start = [mesh.x(9) mesh.y(9) mesh.z(9)];
stop  = [mesh.x(end-8) mesh.y(end-8) mesh.z(end-8)];
[CSX nf2ff] = CreateNF2FFBox(CSX, 'nf2ff', start, stop, [1 1 1 1 0 1]);

%% prepare simulation folder
Sim_Path = 'tmp';
Sim_CSX = 'horn_ant.xml';

[status, message, messageid] = rmdir( Sim_Path, 's' ); % clear previous directory
[status, message, messageid] = mkdir( Sim_Path ); % create empty simulation folder

%% write openEMS compatible xml-file
WriteOpenEMS( [Sim_Path '/' Sim_CSX], FDTD, CSX );

%% show the structure
CSXGeomPlot( [Sim_Path '/' Sim_CSX] );

%% run openEMS
RunOpenEMS( Sim_Path, Sim_CSX);

%% postprocessing & do the plots
U = ReadUI( 'ut1', Sim_Path, freq ); % time domain/freq domain voltage
I = ReadUI( 'it1', Sim_Path, freq ); % time domain/freq domain current (half time step is corrected)

% plot reflection coefficient S11
figure
uf_inc = 0.5*(U.FD{1}.val + I.FD{1}.val .* ZL_a);
if_inc = 0.5*(I.FD{1}.val + U.FD{1}.val ./ ZL_a);
uf_ref = U.FD{1}.val - uf_inc;
if_ref = if_inc - I.FD{1}.val;
s11 = uf_ref ./ uf_inc;
plot( freq/1e9, 20*log10(abs(s11)), 'k-', 'Linewidth', 2 );
ylim([-60 0]);
grid on
title( 'reflection coefficient S_{11}' );
xlabel( 'frequency f / GHz' );
ylabel( 'reflection coefficient |S_{11}|' );

P_in = 0.5*uf_inc .* conj( if_inc ); % antenna feed power

drawnow

%% NFFF contour plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate the far field at phi=0 degrees and at phi=90 degrees
thetaRange = (0:2:359) - 180;
r = 1; % evaluate fields at radius r
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
D_log = 20*log10(nf2ff.E_norm{1}/max(max(nf2ff.E_norm{1})));
% directivity
D_log = D_log + 10*log10(nf2ff.Dmax);

% display polar plot
figure
plot( nf2ff.theta, D_log(:,1) ,'k-' );
xlabel( 'theta (deg)' );
ylabel( 'directivity (dBi)');
grid on;
hold on;
plot( nf2ff.theta, D_log(:,2) ,'r-' );
legend('phi=0','phi=90')

drawnow

%% calculate 3D pattern
phiRange = sort( unique( [-180:5:-100 -100:2.5:-50 -50:1:50 50:2.5:100 100:5:180] ) );
thetaRange = sort( unique([ 0:1:50 50:2.:100 100:5:180 ]));

disp( 'calculating 3D far field...' );
nf2ff = CalcNF2FF(nf2ff, Sim_Path, f0, thetaRange*pi/180, phiRange*pi/180, 'Verbose',2,'Outfile','nf2ff_3D.h5');

E_far_normalized = nf2ff.E_norm{1} / max(nf2ff.E_norm{1}(:)) * nf2ff.Dmax;

[theta,phi] = ndgrid(thetaRange/180*pi,phiRange/180*pi);
x = E_far_normalized .* sin(theta) .* cos(phi);
y = E_far_normalized .* sin(theta) .* sin(phi);
z = E_far_normalized .* cos(theta);
figure
surf( x,y,z, E_far_normalized, 'EdgeColor','none' );
axis equal
axis off
xlabel( 'x' );
ylabel( 'y' );
zlabel( 'z' );

%%
DumpFF2VTK([Sim_Path '/Conical_Horn_Pattern.vtk'],E_far_normalized,thetaRange,phiRange,1e-3);
