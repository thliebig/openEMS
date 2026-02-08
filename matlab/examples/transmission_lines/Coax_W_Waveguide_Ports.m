% Coax_W_Waveguide_Ports.m
% Example demonstrating how to use mode files for waveguide ports, with
% a coaxial T-line.
% (c) 2026 Gadi Lahav (gadi@rfwithcare.com)


clc;
clear;
close all;

function m = snapMesh(m,x0,x1)
    m.x = unique([m.x x0(1) x1(1)]);
    m.x = unique([m.x x0(1) x1(1)]);
    m.y = unique([m.y x0(2) x1(2)]);
    m.y = unique([m.y x0(2) x1(2)]);
    m.z = unique([m.z x0(3) x1(3)]);
    m.z = unique([m.z x0(3) x1(3)]);
end

% -------------------------------------------------------------------------
% Simulation Setup

unit = 1e-3; % drawing unit in mm (adjust if necessary)

% Frequency range and FDTD settings
f_start = 2e9;
f_stop  = 3e9;

% Coax geometry parameters
%coax_D = 2 * unit;                  
%coax_shield_thick = 0.15 * unit;
%coax_wire_D = 0.5 * unit;
%coax_L = 25 * unit;
%
%Airbox_Add = 0 * unit;

coax_D = 2;                  
coax_shield_thick = 0.15;
coax_wire_D = 0.5;
coax_L = 25;

Airbox_Add = 0;


teflon_epsR = 2.5

% -------------------------------------------------------------------------
% Initial simulation definitions

warning off;

addpath('~/opt/openEMS/share/CSXCAD/matlab');
addpath('~/opt/openEMS/share/openEMS/matlab');

FDTD = InitFDTD('NrTS', 30000, 'EndCriteria', 1e-4);
f0 = 0.5*(f_start + f_stop);
FDTD = SetGaussExcite(FDTD, f0, 0.5*(f_stop - f_start));

% Boundary conditions (PML on all sides)
BC = {'PEC','PEC','PEC','PEC','MUR','MUR'};
FDTD = SetBoundaryCond(FDTD, BC);

C0 = 3e8;
mesh_res = ((C0/f_stop)/unit)/100;

% -------------------------------------------------------------------------
% Initial mesh bounding box

mesh.x = [-1 1]*(coax_D*0.5 + coax_shield_thick + Airbox_Add);
mesh.y = [-1 1]*(coax_D*0.5 + coax_shield_thick + Airbox_Add);
mesh.z = [0 coax_L] + [-1 1]*Airbox_Add;

CSX = InitCSX();

% -------------------------------------------------------------------------
% Geometry: Coax & Waveguide port

% PEC conductor
CSX = AddMetal(CSX, 'PEC');

% Coax Outer conductor
CSX = AddCylindricalShell(CSX, 'PEC', 11, [0 0 0], [0 0 coax_L], coax_D*0.5 + coax_shield_thick*0.5, coax_shield_thick);

% Coax Inner conductor
CSX = AddCylinder(CSX, 'PEC', 11, [0 0 0], [0 0 coax_L], coax_wire_D*0.5);

% Teflon Insulator
CSX = AddMaterial( CSX, 'PTFE' );
CSX = SetMaterialProperty( CSX, 'PTFE', 'Epsilon', teflon_epsR );

% Coax Inner conductor
CSX = AddCylinder(CSX, 'PTFE', 10, [0 0 0], [0 0 coax_L], coax_D*0.5);

% Add all relevant points to mesh
mesh.x = [mesh.x linspace(-coax_D*0.5,-coax_wire_D*0.5,7)];
mesh.x = [mesh.x linspace(-(coax_D*0.5 + coax_shield_thick),-coax_D*0.5,4)];
mesh.x = [mesh.x linspace(-coax_wire_D*0.5,coax_wire_D*0.5,6)];
mesh.x = [mesh.x linspace(coax_D*0.5,(coax_D*0.5 + coax_shield_thick),4)];
mesh.x = [mesh.x linspace(coax_wire_D*0.5,coax_D*0.5,7)];

mesh.x = unique(mesh.x);
mesh.y = mesh.x;
mesh.z = [0 coax_L];

% Smooth all meshes
mesh.x = SmoothMeshLines(mesh.x, mesh_res, 1.25);
mesh.y = SmoothMeshLines(mesh.y, mesh_res, 1.25);
mesh.z = SmoothMeshLines(mesh.z, mesh_res, 1.25);

% -------------------------------------------------------------------------
% Ports Setup
%
%% --- Coaxial Port ---
%% Coax start and stop for lumped port
%start_coax = [0 0 0];
%stop_coax  = [0 0 coax_length];
%
%% Add coaxial lumped port
%[CSX, port{1}] = AddLumpedPort(CSX, 5, 1, 50, start_coax, stop_coax, [0 0 1], true);
%
%% --- Rectangular Waveguide Port ---
%% Note: If AddWaveguidePort isn't fully implemented, leave placeholder
%start_wg = [-wg_a/2 -wg_b/2 coax_length];
%stop_wg  = [ wg_a/2  wg_b/2 coax_length+1e-6]; % small extension
%
%% Placeholder for your AddWaveguidePort version
%% [CSX, port{2}] = AddWaveguidePort(CSX, 5, 2, start_wg, stop_wg, ...
%%                                   dir, E_func, H_func, kc, excite);
%% For now, you can use rectangular port if available:
%mode_name = 'TE10'; % example rectangular waveguide mode
%[CSX, port{2}] = AddRectWaveGuidePort(CSX, 5, 2, start_wg, stop_wg, 'z', ...
%                                      wg_a, wg_b, mode_name, 1);

% -------------------------------------------------------------------------
% Write XML and Run Simulation

% Define the grid inside CSX
CSX = DefineRectGrid(CSX, unit, mesh);

Sim_Path = 'tmp';
Sim_CSX  = 'coax_wg.xml';

WriteOpenEMS([Sim_Path '/' Sim_CSX], FDTD, CSX);

% show the structure
CSXGeomPlot( [Sim_Path '/' Sim_CSX] );

%!openEMS [Sim_Path '/' Sim_CSX];

% -------------------------------------------------------------------------
% Post-processing

freq = linspace(f_start, f_stop, 201);

port = calcPort(port, Sim_Path, freq);

% Calculate and plot S-parameters
S11 = port{1}.uf.ref ./ port{1}.uf.inc;
S21 = port{2}.uf.ref ./ port{1}.uf.inc;

figure;
plot(freq/1e9, 20*log10(abs(S11)), 'b-', 'LineWidth', 2);
hold on;
plot(freq/1e9, 20*log10(abs(S21)), 'r--', 'LineWidth', 2);
xlabel('Frequency (GHz)');
ylabel('S-Parameters (dB)');
legend('S11', 'S21');
title('Coax to Waveguide Test (Matlab)');
grid on;


