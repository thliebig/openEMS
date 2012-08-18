function directional_coupler
%
% EXAMPLE / microstrip / directional_coupler
%
% Stacked directional coupler in microstrip technology.
%
% This example demonstrates:
%  - simple microstrip geometry
%  - S-parameter calculation using the ypar-method
%  - display of coupler parameters
%  - display of S11 (smith chart)
% 
%
% Tested with
%  - Matlab 2010b
%  - Octave 3.2.4
%  - openEMS v0.0.17
%
% (C) 2010 Sebastian Held <sebastian.held@gmx.de>

clear
close all
clc

% sim settings
showStructure = 1;
runSimulation = 1;

for n=1:4
    if n > 1, showStructure = 0; end
    ports{n} = sim( n, showStructure, runSimulation );
end
postprocess( ports );




function ports = sim( simnr, showStructure, runSimulation )
physical_constants

% setup the simulation
drawingunit = 1e-6; % specify everything in um
Sim_Path = ['tmp' int2str(simnr)];
Sim_CSX  = 'tmp.xml';
f_max    = 100e6;
lambda   = c0/f_max;

% specify the coupler
pcb1.w   = 147000;
pcb1.h   = 54500;
pcb1.t   = 1524;
pcb1.epr = 3;
msl1.w   = 135000;
msl1.h   = 2800;
pcb2.w   = 107000;
pcb2.h   = 14000;
pcb2.t   = 1524;
pcb2.epr = 3;
msl2.w   = 95000;
msl2.h   = 4000;


CSX = InitCSX();

% create the mesh
mesh.x = [-pcb1.w/2 pcb1.w/2 -pcb2.w/2 pcb2.w/2 -msl1.w/2 msl1.w/2 -msl2.w/2 msl2.w/2];
mesh.x = [mesh.x linspace(-msl2.w/2,-msl2.w/2+msl2.h, 5) linspace(msl2.w/2,msl2.w/2-msl2.h, 5)];
mesh.y = [-pcb1.h/2 pcb1.h/2 -pcb2.h/2 pcb2.h/2 -msl1.h/2 msl1.h/2 -msl2.h/2 msl2.h/2];
mesh.z = [linspace(0,pcb1.t,5) linspace(pcb1.t,pcb1.t+pcb2.t,5)];
mesh.z = [mesh.z mesh.z(end)+10*(mesh.z(end)-mesh.z(1))]; % add space above pcb
res = lambda/sqrt(max([pcb1.epr,pcb2.epr])) / 20 / drawingunit;
mesh.x = SmoothMeshLines2(mesh.x,res);
mesh.y = SmoothMeshLines2(mesh.y,res);
mesh.z = SmoothMeshLines2(mesh.z,res);
mesh = AddPML( mesh, [8 8 8 8 8 8] ); % add space for PML
CSX = DefineRectGrid( CSX, drawingunit, mesh );

%% create the structure

% microstrip
CSX = AddMetal( CSX, 'PEC' );
start = [-msl1.w/2, -msl1.h/2, pcb1.t];
stop  = [ msl1.w/2,  msl1.h/2, pcb1.t];
priority = 100; % the geometric priority is set to 100
CSX = AddBox( CSX, 'PEC', priority, start, stop );

% ground plane
CSX = AddMetal( CSX, 'PEC_ground' );
start = [-pcb1.w/2, -pcb1.h/2, 0];
stop  = [ pcb1.w/2,  pcb1.h/2, 0];
CSX = AddBox( CSX, 'PEC_ground', priority, start, stop );

% substrate 1
start = [-pcb1.w/2, -pcb1.h/2, 0];
stop  = [ pcb1.w/2,  pcb1.h/2, pcb1.t];
priority = 10;
CSX = AddMaterial( CSX, 'substrate1' );
CSX = SetMaterialProperty( CSX, 'substrate1', 'Epsilon', pcb1.epr );
CSX = AddBox( CSX, 'substrate1', priority, start, stop );

% substrate 2
start = [-pcb2.w/2, -pcb2.h/2, pcb1.t];
stop  = [ pcb2.w/2,  pcb2.h/2, pcb1.t+pcb2.t];
priority = 10;
CSX = AddMaterial( CSX, 'substrate2' );
CSX = SetMaterialProperty( CSX, 'substrate2', 'Epsilon', pcb2.epr );
CSX = AddBox( CSX, 'substrate2', priority, start, stop );

% stripline
start = [-msl2.w/2, -msl2.h/2, pcb1.t+pcb2.t];
stop  = [ msl2.w/2,  msl2.h/2, pcb1.t+pcb2.t];
priority = 100;
CSX = AddBox( CSX, 'PEC', priority, start, stop );

% connections
start = [-msl2.w/2,        -msl2.h/2, pcb1.t+pcb2.t];
stop  = [-msl2.w/2+msl2.h, -pcb2.h/2, pcb1.t+pcb2.t];
priority = 100;
CSX = AddBox( CSX, 'PEC', priority, start, stop );
start = [ msl2.w/2,        -msl2.h/2, pcb1.t+pcb2.t];
stop  = [ msl2.w/2-msl2.h, -pcb2.h/2, pcb1.t+pcb2.t];
priority = 100;
CSX = AddBox( CSX, 'PEC', priority, start, stop );

%% ports
% this project needs 4 simulations
for n=1:4
    portexcite{n} = [];
end
portexcite{simnr} = 'excite';

% port 1: input port
start = [-msl1.w/2, 0, pcb1.t];
stop  = [-msl1.w/2, 0, 0];
[CSX ports{1}] = AddCurvePort( CSX, 999, 1, 50, start, stop, portexcite{1} );
% port 2: output port
start = [msl1.w/2, 0, pcb1.t];
stop  = [msl1.w/2, 0, 0];
[CSX ports{2}] = AddCurvePort( CSX, 999, 2, 50, start, stop, portexcite{2} );
% port 3: coupled port
start = [-msl2.w/2+msl2.h/2, -pcb2.h/2, pcb1.t+pcb2.t];
stop  = [-msl2.w/2+msl2.h/2, -pcb2.h/2, 0];
[CSX ports{3}] = AddCurvePort( CSX, 999, 3, 50, start, stop, portexcite{3} );
% port 4: isolated port
start = [msl2.w/2-msl2.h/2, -pcb2.h/2, pcb1.t+pcb2.t];
stop  = [msl2.w/2-msl2.h/2, -pcb2.h/2, 0];
[CSX ports{4}] = AddCurvePort( CSX, 999, 4, 50, start, stop, portexcite{4} );

%% setup FDTD parameters & excitation function
max_timesteps = 50000;
min_decrement = 1e-6;
FDTD = InitFDTD( max_timesteps, min_decrement );
FDTD = SetGaussExcite( FDTD, 0, f_max );
BC = {'PML_8' 'PML_8' 'PML_8' 'PML_8' 'PML_8' 'PML_8'};
BC = {'MUR' 'MUR' 'MUR' 'MUR' 'MUR' 'MUR'}; % faster
FDTD = SetBoundaryCond( FDTD, BC );

%% Write openEMS compatible xml-file
if runSimulation
    [dummy,dummy,dummy] = rmdir(Sim_Path,'s');
end
[dummy,dummy,dummy] = mkdir(Sim_Path);
WriteOpenEMS([Sim_Path '/' Sim_CSX],FDTD,CSX);

if showStructure
    CSXGeomPlot( [Sim_Path '/' Sim_CSX] );
end

%% run openEMS
openEMS_opts = '';
openEMS_opts = [openEMS_opts ' --engine=fastest'];
% openEMS_opts = [openEMS_opts ' --debug-material'];
% openEMS_opts = [openEMS_opts ' --debug-boxes'];
if runSimulation
    RunOpenEMS( Sim_Path, Sim_CSX, openEMS_opts );
end




function postprocess( ports )
f = linspace( 0, 100e6, 201 );
Y = calc_ypar( f, ports{1}, 'tmp' );
R = 50;
S = y2s(Y,R);

% insertion loss
IL_dB = -20 * log10(abs(squeeze(S(2,1,:))));

% coupling factor
CF_dB = -20 * log10(abs(squeeze(S(3,1,:))));

% isolation
I_dB  = -20 * log10(abs(squeeze(S(4,1,:))));

% directivity
D_dB  = -20 * log10(abs(squeeze(S(4,1,:) ./ S(3,1,:))));

figure
plot( f, [IL_dB CF_dB I_dB D_dB] );
legend( {'insertion loss','coupling factor','isolation','directivity'} );
title( ['performance of the coupler for a termination resistance of R=' num2str(R)] );
grid on

smithchart
S11 = squeeze(S(1,1,:));
plot( real(S11), imag(S11) );
legend( 'S_{11}' );
title( ['performance of the coupler for a termination resistance of R=' num2str(R)] );
axis( [-1 1 -1 1] );



function smithchart
% smith chart
figure
if exist( 'smith', 'file' )
    % smith chart
    % www.ece.rutgers.edu/~orfanidi/ewa
    % or cmt toolbox from git.ate.uni-duisburg.de
    smith
else
    % poor man smith chart
    color = [.6 .6 .6];
    h = plot( sin(0:0.01:2*pi), cos(0:0.01:2*pi), 'Color', color );
    hg = hggroup;
    set( h,'Parent',hg );
    hold on
    plot( hg, 0.25+0.75*sin(0:0.01:2*pi), 0.75*cos(0:0.01:2*pi), 'Color', color );
    plot( hg, 0.5+0.5*sin(0:0.01:2*pi), 0.5*cos(0:0.01:2*pi), 'Color', color );
    plot( hg, 0.75+0.25*sin(0:0.01:2*pi), 0.25*cos(0:0.01:2*pi), 'Color', color );
    plot( hg, [-1 1], [0 0], 'Color', color );
    axis equal
    axis off
end


function s = y2s(y, ZL)
% S = y2s(Y, ZL)
%
% Admittance to Scattering transformation
% for square matrices at multiple frequencies
%
% ZL defaults to 50 Ohm

if nargin < 2
    ZL = 50;
end

if size(size(y),2) > 2
    nF = size(y,3);
else
    nF = 1;
end

I = diag(ones(1, size(y,2)))/ZL;

for i=1:nF
    %s(:,:,i) = inv(I+y(:,:,i)) * (I-y(:,:,i));
    s(:,:,i) = (I+y(:,:,i)) \ (I-y(:,:,i));
end
