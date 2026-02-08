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

Mode_File_Path = '../../../python/Tests';
E_Mode_File = 'Coax_Er.csv';
H_Mode_File = 'Coax_Hr.csv';

teflon_epsR = 2.5;

% -------------------------------------------------------------------------
% Initial simulation definitions

warning off;

addpath('~/opt/openEMS/share/CSXCAD/matlab');
addpath('~/opt/openEMS/share/openEMS/matlab');

physical_constants;

FDTD = InitFDTD('NrTS', 30000, 'EndCriteria', 1e-4);
f0 = 0.5*(f_start + f_stop);
FDTD = SetGaussExcite(FDTD, f0, 0.5*(f_stop - f_start));

% Boundary conditions (PML on all sides)
BC = {'PEC','PEC','PEC','PEC','MUR','MUR'};
FDTD = SetBoundaryCond(FDTD, BC);

mesh_res = ((c0/f_stop)/unit)/100;

% -------------------------------------------------------------------------
% Initial mesh bounding box


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
mesh.x = [];
mesh.y = [];
mesh.z = [];
mesh.x = [mesh.x linspace(-coax_D*0.5,-coax_wire_D*0.5,7)];
mesh.x = [mesh.x linspace(-(coax_D*0.5 + coax_shield_thick),-coax_D*0.5,4)];
mesh.x = [mesh.x linspace(-coax_wire_D*0.5,coax_wire_D*0.5,6)];
mesh.x = [mesh.x linspace(coax_D*0.5,(coax_D*0.5 + coax_shield_thick),4)];
mesh.x = [mesh.x linspace(coax_wire_D*0.5,coax_D*0.5,7)];

mesh.x = unique(mesh.x);
mesh.y = mesh.x;
mesh.z = unique([0:mesh_res*0.25:mesh_res (coax_L - fliplr(0:mesh_res*0.25:mesh_res)) coax_L]);



% Smooth all meshes
mesh.x = SmoothMeshLines(mesh.x, mesh_res, 1.25);
mesh.y = SmoothMeshLines(mesh.y, mesh_res, 1.25);
mesh.z = SmoothMeshLines(mesh.z, mesh_res, 1.25);


% Define the grid inside CSX
CSX = DefineRectGrid(CSX, unit, mesh);

% -------------------------------------------------------------------------
% Ports Setup

% --- Coaxial Port ---

% Wave number
kz = 82.84554871;
% Wave impedance
Zw = 238.26517157;
% Line impedance
Zl = 52.43928218;

% Calculate cutoff wave number 
kc = sqrt((2*pi*f0/c0)^2 - kz^2);

% Add coaxial lumped port
% Coax start and stop for lumped port
start_coax = [(coax_D*0.5 + coax_shield_thick)*[-1 -1] mesh.z(2)];
stop_coax  = [(coax_D*0.5 + coax_shield_thick)*[ 1  1] mesh.z(3)];

[CSX, port{1}] = AddWaveGuidePort(CSX, 15, 1, start_coax, stop_coax, 'z', '', '', kc , 1, 'E_WG_file', E_Mode_File, 'H_WG_file', H_Mode_File);

% Add coaxial lumped port
% Coax start and stop for lumped port
start_coax = [(coax_D*0.5 + coax_shield_thick)*[-1 -1] mesh.z(end - 1)];
stop_coax  = [(coax_D*0.5 + coax_shield_thick)*[ 1  1] mesh.z(end - 2)];

[CSX, port{2}] = AddWaveGuidePort(CSX, 15, 2, start_coax, stop_coax, 'z', '', '', kc , 1, 'E_WG_file', E_Mode_File, 'H_WG_file', H_Mode_File);
% -------------------------------------------------------------------------
% Write XML and Run Simulation


Sim_Path = 'tmp';
Sim_CSX  = 'coax_wg.xml';

WriteOpenEMS([Sim_Path '/' Sim_CSX], FDTD, CSX);

copyfile([Mode_File_Path '/' E_Mode_File], [Sim_Path '/']);
copyfile([Mode_File_Path '/' H_Mode_File], [Sim_Path '/']);

% show the structure
CSXGeomPlot( [Sim_Path '/' Sim_CSX] );

% run openEMS
RunOpenEMS( Sim_Path, Sim_CSX );

% -------------------------------------------------------------------------
% Post-processing

freq = linspace(f_start, f_stop, 201);

port = calcPort(port, Sim_Path, freq, 'RefImpedance', 50);

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


