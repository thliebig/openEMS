%
% EXAMPLE / antennas / inverted-f antenna (ifa) 2.4GHz
%
% This example demonstrates how to:
%  - calculate the reflection coefficient of an ifa
%  - calculate farfield of an ifa
%
% Tested with
%  - Octave 3.7.5
%  - openEMS v0.0.30+ (git 10.07.2013)
%
% (C) 2013 Stefan Mahr <dac922@gmx.de>

close all
clear
clc

%% setup the simulation
physical_constants;
unit = 1e-3; % all length in mm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                substrate.width
%  _______________________________________________    __ substrate.
% | A                        ifa.l                |\  __    thickness
% | |ifa.e         __________________________     | |
% | |             |    ___  _________________| w2 | |
% | |       ifa.h |   |   ||                      | |
% |_V_____________|___|___||______________________| |
% |                .w1   .wf\                     | |
% |                   |.fp|  \                    | |
% |                       |    feed point         | |
% |                       |                       | | substrate.length
% |<- substrate.width/2 ->|                       | |
% |                                               | |
% |_______________________________________________| |
%  \_______________________________________________\|
%
% Note: It's not checked whether your settings make sense, so check
%       graphical output carefully.
%
substrate.width  = 80;             % width of substrate
substrate.length = 80;             % length of substrate
substrate.thickness = 1.5;         % thickness of substrate
substrate.cells = 4;               % use 4 cells for meshing substrate

ifa.h  = 8;            % height of short circuit stub
ifa.l  = 22.5;         % length of radiating element
ifa.w1 = 4;            % width of short circuit stub
ifa.w2 = 2.5;          % width of radiating element
ifa.wf = 1;            % width of feed element
ifa.fp = 4;            % position of feed element relative to short
                       %  circuit stub
ifa.e  = 10;           % distance to edge


% substrate setup
substrate.epsR   = 4.3;
substrate.kappa  = 1e-3 * 2*pi*2.45e9 * EPS0*substrate.epsR;

%setup feeding
feed.R = 50;     %feed resistance

%open AppCSXCAD and show ifa
show = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% size of the simulation box
SimBox = [substrate.width*2 substrate.length*2 150];

%% setup FDTD parameter & excitation function
f0 = 2.5e9; % center frequency
fc = 1e9; % 20 dB corner frequency

FDTD = InitFDTD('NrTS',  60000 );
FDTD = SetGaussExcite( FDTD, f0, fc );
BC = {'MUR' 'MUR' 'MUR' 'MUR' 'MUR' 'MUR'}; % boundary conditions
FDTD = SetBoundaryCond( FDTD, BC );

%% setup CSXCAD geometry & mesh
CSX = InitCSX();

%initialize the mesh with the "air-box" dimensions
mesh.x = [-SimBox(1)/2 SimBox(1)/2];
mesh.y = [-SimBox(2)/2 SimBox(2)/2];
mesh.z = [-SimBox(3)/2 SimBox(3)/2];

%% create substrate
CSX = AddMaterial( CSX, 'substrate');
CSX = SetMaterialProperty( CSX, 'substrate', 'Epsilon',substrate.epsR, 'Kappa', substrate.kappa);
start = [-substrate.width/2  -substrate.length/2                    0];
stop  = [ substrate.width/2   substrate.length/2  substrate.thickness];
CSX = AddBox( CSX, 'substrate', 1, start, stop );
% add extra cells to discretize the substrate thickness
mesh.z = [linspace(0,substrate.thickness,substrate.cells+1) mesh.z];

%% create ground plane
CSX = AddMetal( CSX, 'groundplane' ); % create a perfect electric conductor (PEC)
start = [-substrate.width/2  -substrate.length/2        substrate.thickness];
stop  = [ substrate.width/2   substrate.length/2-ifa.e  substrate.thickness];
CSX = AddBox(CSX, 'groundplane', 10, start,stop);

%% create ifa
CSX = AddMetal( CSX, 'ifa' ); % create a perfect electric conductor (PEC)
tl = [0,substrate.length/2-ifa.e,substrate.thickness];   % translate
start = [0 0.5 0] + tl;
stop = start + [ifa.wf ifa.h-0.5 0];
CSX = AddBox( CSX, 'ifa', 10,  start, stop);  % feed element
start = [-ifa.fp 0 0] + tl;
stop =  start + [-ifa.w1 ifa.h 0];
CSX = AddBox( CSX, 'ifa', 10,  start, stop);  % short circuit stub
start = [(-ifa.fp-ifa.w1) ifa.h 0] + tl;
stop = start + [ifa.l -ifa.w2 0];
CSX = AddBox( CSX, 'ifa', 10, start, stop);   % radiating element

ifa_mesh = DetectEdges(CSX, [], 'SetProperty','ifa');
mesh.x = [mesh.x SmoothMeshLines(ifa_mesh.x, 0.5)];
mesh.y = [mesh.y SmoothMeshLines(ifa_mesh.y, 0.5)];

%% apply the excitation & resist as a current source
start = [0 0 0] + tl;
stop  = start + [ifa.wf 0.5 0];
[CSX port] = AddLumpedPort(CSX, 5 ,1 ,feed.R, start, stop, [0 1 0], true);

%% finalize the mesh
% generate a smooth mesh with max. cell size: lambda_min / 20
mesh = DetectEdges(CSX, mesh);
mesh = SmoothMesh(mesh, c0 / (f0+fc) / unit / 20);
CSX = DefineRectGrid(CSX, unit, mesh);

%% add a nf2ff calc box; size is 3 cells away from MUR boundary condition
start = [mesh.x(4)     mesh.y(4)     mesh.z(4)];
stop  = [mesh.x(end-3) mesh.y(end-3) mesh.z(end-3)];
[CSX nf2ff] = CreateNF2FFBox(CSX, 'nf2ff', start, stop);

%% prepare simulation folder
Sim_Path = 'tmp_IFA';
Sim_CSX = 'IFA.xml';

try confirm_recursive_rmdir(false,'local'); end
 
[status, message, messageid] = rmdir( Sim_Path, 's' ); % clear previous directory
[status, message, messageid] = mkdir( Sim_Path ); % create empty simulation folder

%% write openEMS compatible xml-file
WriteOpenEMS( [Sim_Path '/' Sim_CSX], FDTD, CSX );

%% show the structure
if (show == 1)
  CSXGeomPlot( [Sim_Path '/' Sim_CSX] );
end


%% run openEMS
RunOpenEMS( Sim_Path, Sim_CSX);  %RunOpenEMS( Sim_Path, Sim_CSX, '--debug-PEC -v');

%% postprocessing & do the plots
freq = linspace( max([1e9,f0-fc]), f0+fc, 501 );
port = calcPort(port, Sim_Path, freq);

Zin = port.uf.tot ./ port.if.tot;
s11 = port.uf.ref ./ port.uf.inc;
P_in = real(0.5 * port.uf.tot .* conj( port.if.tot )); % antenna feed power

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
plot( freq/1e6, 20*log10(abs(s11)), 'k-', 'Linewidth', 2 );
grid on
title( 'reflection coefficient S_{11}' );
xlabel( 'frequency f / MHz' );
ylabel( 'reflection coefficient |S_{11}|' );

drawnow

%% NFFF contour plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find resonance frequency from s11
f_res_ind = find(s11==min(s11));
f_res = freq(f_res_ind);

%%
disp( 'calculating 3D far field pattern and dumping to vtk (use Paraview to visualize)...' );
thetaRange = (0:2:180);
phiRange = (0:2:360) - 180;
nf2ff = CalcNF2FF(nf2ff, Sim_Path, f_res, thetaRange*pi/180, phiRange*pi/180,'Verbose',1,'Outfile','3D_Pattern.h5');

plotFF3D(nf2ff)

% display power and directivity
disp( ['radiated power: Prad = ' num2str(nf2ff.Prad) ' Watt']);
disp( ['directivity: Dmax = ' num2str(nf2ff.Dmax) ' (' num2str(10*log10(nf2ff.Dmax)) ' dBi)'] );
disp( ['efficiency: nu_rad = ' num2str(100*nf2ff.Prad./real(P_in(f_res_ind))) ' %']);

E_far_normalized = nf2ff.E_norm{1} / max(nf2ff.E_norm{1}(:)) * nf2ff.Dmax;
DumpFF2VTK([Sim_Path '/3D_Pattern.vtk'],E_far_normalized,thetaRange,phiRange,1e-3);
