function [port nf2ff] = Patch_Antenna_Array(Sim_Path, postproc_only, show_structure, xpos, caps, resist, active )
% [port nf2ff] = Patch_Antenna_Array(Sim_Path, postproc_only, show_structure, xpos, caps, resist, active )
%
% Script to setup the patch array as described in [1].
% Run main script in Patch_Antenna_Phased_Array.m instead!
%
% Sim_Path: Simulation path
% postproc_only: set to post process only 0/1
% show_structure: show the structure in AppCSXCAD 0/1
% xpos: the x-position for each antenna is defined
% caps: the port capacity (will override active port)
% resist: port resistance
% active: switch port active
%
% References:
% [1] Y. Yusuf and X. Gong, “A low-cost patch antenna phased array with
%   analog beam steering using mutual coupling and reactive loading,” IEEE
%   Antennas Wireless Propag. Lett., vol. 7, pp. 81–84, 2008.
%
% Tested with
%  - Matlab 2011a
%  - openEMS v0.0.31
%
% (C) 2013 Thorsten Liebig <thorsten.liebig@gmx.de>

% example
% xpos = [-41 0 41];
% caps = [0.2e-12 0 0.2e-12];
% active = [0 1 0];
% resist = [50 50 50];

%% setup the simulation
physical_constants;
unit = 1e-3; % all length in mm

% patch geometry setup
patch.W  = 35;  % width
patch.L = 28.3; % length
patch.Ws = 3.8; % width of feeding stub
patch.Gs = 1;   % width of feeding gab
patch.l = 6;    % length of feeding stub
patch.y0 = 10;  % depth of feeding stub into into patch

% patch resonance frequency
f0 = 3e9;

%substrate setup
substrate.name = 'Ro3003';
substrate.epsR   = 3;
substrate.kappa  = 0.0013 * 2*pi*f0 * EPS0*substrate.epsR;
substrate.thickness = 1.524;
substrate.cells = 4;

substrate.width = patch.W + max(xpos) - min(xpos) + 4*patch.l;
substrate.length = 3*patch.l + patch.L;

% size of the simulation box
AirSpacer = [50 50 30];

edge_res = [-1/3 2/3]*1;

%% setup FDTD parameter & excitation function
fc = 2e9; % 20 dB corner frequency
FDTD = InitFDTD( 'EndCriteria', 1e-4 );
FDTD = SetGaussExcite( FDTD, f0, fc );
BC = [1 1 1 1 1 1]*3;
FDTD = SetBoundaryCond( FDTD, BC );

%% setup CSXCAD geometry & mesh
CSX = InitCSX();

mesh.x = [];
mesh.y = [];
mesh.z = [];

%% create patch
CSX = AddMetal( CSX, 'patch' ); % create a perfect electric conductor (PEC)

for port_nr=1:numel(xpos)
    start = [xpos(port_nr)-patch.W/2           patch.l         substrate.thickness];
    stop  = [xpos(port_nr)-patch.Ws/2-patch.Gs patch.l+patch.L substrate.thickness];
    CSX = AddBox(CSX,'patch',10, start, stop);
    mesh.x = [mesh.x xpos(port_nr)-patch.W/2-edge_res];

    start = [xpos(port_nr)+patch.W/2           patch.l         substrate.thickness];
    stop  = [xpos(port_nr)+patch.Ws/2+patch.Gs patch.l+patch.L substrate.thickness];
    CSX = AddBox(CSX,'patch',10, start, stop);
    mesh.x = [mesh.x xpos(port_nr)+patch.W/2+edge_res];

    mesh.y = [mesh.y patch.l-edge_res patch.l+patch.L+edge_res];

    start = [xpos(port_nr)-patch.Ws/2-patch.Gs patch.l+patch.y0 substrate.thickness];
    stop  = [xpos(port_nr)+patch.Ws/2+patch.Gs patch.l+patch.L  substrate.thickness];
    CSX = AddBox(CSX,'patch',10, start, stop);

    % feed line
    start = [xpos(port_nr)-patch.Ws/2 patch.l+patch.y0 substrate.thickness];
    stop  = [xpos(port_nr)+patch.Ws/2 0                substrate.thickness];
    CSX = AddBox(CSX,'patch',10, start, stop);

    mesh.x = [mesh.x xpos(port_nr)+linspace(-patch.Ws/2-patch.Gs,-patch.Ws/2,3) xpos(port_nr)+linspace(patch.Ws/2,patch.Ws/2+patch.Gs,3)];

    start = [xpos(port_nr)-patch.Ws/2 0 0];
    stop  = [xpos(port_nr)+patch.Ws/2 0 substrate.thickness];
    if (caps(port_nr)>0)
        CSX = AddLumpedElement(CSX, ['C_' num2str(port_nr)], 2, 'C', caps(port_nr));
        CSX = AddBox(CSX,['C_' num2str(port_nr)],10, start, stop);

        [CSX port{port_nr}] = AddLumpedPort(CSX, 5 ,port_nr ,inf, start, stop, [0 0 1], 0);
    else
        % feed port
        [CSX port{port_nr}] = AddLumpedPort(CSX, 5 ,port_nr, resist(port_nr), start, stop, [0 0 1], active(port_nr));
    end
end

%% create substrate
CSX = AddMaterial( CSX, substrate.name );
CSX = SetMaterialProperty( CSX, substrate.name, 'Epsilon', substrate.epsR, 'Kappa', substrate.kappa );
start = [-substrate.width/2 0                0];
stop  = [ substrate.width/2 substrate.length substrate.thickness];
CSX = AddBox( CSX, substrate.name, 0, start, stop );

mesh.x = [mesh.x start(1) stop(1)];
mesh.y = [mesh.y start(2) stop(2)];

% add extra cells to discretize the substrate thickness
mesh.z = [linspace(0,substrate.thickness,substrate.cells+1) mesh.z];

%% create ground (same size as substrate)
CSX = AddMetal( CSX, 'gnd' ); % create a perfect electric conductor (PEC)
start(3)=0;
stop(3) =0;
CSX = AddBox(CSX,'gnd',10,start,stop);

%% finalize the mesh
% generate a smooth mesh with max. cell size: lambda_min / 20
mesh = SmoothMesh(mesh, 2, 1.3);
mesh.x = [mesh.x min(mesh.x)-AirSpacer(1) max(mesh.x)+AirSpacer(1)];
mesh.y = [mesh.y min(mesh.y)-AirSpacer(2) max(mesh.y)+AirSpacer(2)];
mesh.z = [mesh.z min(mesh.z)-AirSpacer(3) max(mesh.z)+2*AirSpacer(3)];

mesh = SmoothMesh(mesh, c0 / (f0+fc) / unit / 20, 1.3);

%% add a nf2ff calc box; size is 3 cells away from MUR boundary condition
start = [mesh.x(4)     mesh.y(4)     mesh.z(4)];
stop  = [mesh.x(end-3) mesh.y(end-3) mesh.z(end-3)];
[CSX nf2ff] = CreateNF2FFBox(CSX, 'nf2ff', start, stop);

mesh = AddPML(mesh,(BC==3)*8);
CSX = DefineRectGrid(CSX, unit, mesh);

%% prepare simulation folder
Sim_CSX = 'patch_array.xml';

if (postproc_only==0)
    [status, message, messageid] = rmdir( Sim_Path, 's' ); % clear previous directory
    [status, message, messageid] = mkdir( Sim_Path ); % create empty simulation folder

    %% write openEMS compatible xml-file
    WriteOpenEMS( [Sim_Path '/' Sim_CSX], FDTD, CSX );

    %% show the structure
    if (show_structure>0)
        CSXGeomPlot( [Sim_Path '/' Sim_CSX] );
    end

    %% run openEMS
    RunOpenEMS( Sim_Path, Sim_CSX);
end

