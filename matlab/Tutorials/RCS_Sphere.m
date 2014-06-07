%
% Tutorials / radar cross section of a metal sphere
%
% Describtion at:
% http://openems.de/index.php/Tutorial:_RCS_Sphere
%
% Tested with
%  - Matlab 2013a / Octave 3.8.1
%  - openEMS v0.0.32
%
% (C) 2012-2014 Thorsten Liebig <thorsten.liebig@uni-due.de>

close all
clear
clc

%% setup the simulation
physical_constants;
unit = 1e-3; % all length in mm

sphere.rad = 200;

inc_angle = 0 /180*pi; %incident angle (to x-axis) in rad

% size of the simulation box
SimBox = 1000;
PW_Box = 750;

%% setup FDTD parameter & excitation function
f_start =  50e6; % start frequency
f_stop = 1000e6; % stop  frequency
f0 = 500e6;

FDTD = InitFDTD( );
FDTD = SetGaussExcite( FDTD, 0.5*(f_start+f_stop), 0.5*(f_stop-f_start) );
BC = [1 1 1 1 1 1]*3;  % set boundary conditions
FDTD = SetBoundaryCond( FDTD, BC );

%% setup CSXCAD geometry & mesh
max_res = c0 / f_stop / unit / 20; % cell size: lambda/20
CSX = InitCSX();

%create mesh
smooth_mesh = SmoothMeshLines([0 SimBox/2], max_res);
mesh.x = unique([-smooth_mesh smooth_mesh]);
mesh.y = mesh.x;
mesh.z = mesh.x;

%% create metal sphere
CSX = AddMetal( CSX, 'sphere' ); % create a perfect electric conductor (PEC)
CSX = AddSphere(CSX,'sphere',10,[0 0 0],sphere.rad);

%% plane wave excitation
k_dir = [cos(inc_angle) sin(inc_angle) 0]; % plane wave direction
E_dir = [0 0 1]; % plane wave polarization --> E_z

CSX = AddPlaneWaveExcite(CSX, 'plane_wave', k_dir, E_dir, f0);
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
EF = ReadUI( 'et', Sim_Path, f0 ); % time domain/freq domain voltage
Pin = 0.5*norm(E_dir)^2/Z0 .* abs(EF.FD{1}.val).^2;

%%
nf2ff = CalcNF2FF(nf2ff, Sim_Path, f0, pi/2, [-180:2:180]*pi/180, 'Mode',1);
RCS = 4*pi./Pin(1).*nf2ff.P_rad{1}(:);
polar(nf2ff.phi,RCS);
xlabel('x -->');
ylabel('y -->');
hold on
grid on

drawnow

%%
freq = linspace(f_start,f_stop,100);
EF = ReadUI( 'et', Sim_Path, freq ); % time domain/freq domain voltage
Pin = 0.5*norm(E_dir)^2/Z0 .* abs(EF.FD{1}.val).^2;

nf2ff = CalcNF2FF(nf2ff, Sim_Path, freq, pi/2, pi+inc_angle, 'Mode',1);
for fn=1:numel(freq)
    back_scat(fn) = 4*pi./Pin(fn).*nf2ff.P_rad{fn}(1);
end

%%
figure
plot(freq/1e6,back_scat,'Linewidth',2);
grid on;
xlabel('frequency (MHz) \rightarrow');
ylabel('RCS (m^2) \rightarrow');
title('radar cross section');

%%
figure
lambda = c0./freq;
semilogy(sphere.rad*unit./lambda,back_scat/(pi*sphere.rad*unit*sphere.rad*unit),'Linewidth',2);
ylim([10^-2 10^1])
grid on;
xlabel('sphere radius / wavelength \rightarrow');
ylabel('RCS / (\pi a^2) \rightarrow');
title('normalized radar cross section');
