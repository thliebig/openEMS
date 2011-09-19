%
% EXAMPLE / microstrip / MSL2
%
% This example shows how to use the MSL-port.
% The MSL is excited at the center of the computational volume. The
% boundary at xmin is an absorbing boundary (Mur) and at xmax an electric
% wall. The reflection coefficient at this wall is S11 = -1.
% Direction of propagation is x.
%
% This example demonstrates:
%  - simple microstrip geometry (made of PEC)
%  - MSL port
%  - MSL analysis
%
% You may modify the PEC boundary condition at xmax to become a MUR
% boundary. This resembles a matched microstrip line.
%
% Tested with
%  - Matlab 2009b
%  - Octave 3.3.52
%  - openEMS v0.0.14
%
% (C) 2010 Sebastian Held <sebastian.held@uni-due.de>

close all
clear
clc

%% switches
postproc_only = 0;

%% setup the simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
physical_constants;
unit = 1e-6; % specify everything in um
MSL_length    = 10000;
MSL_width     = 1000;
substrate_thickness = 254;
substrate_epr = 3.66;

% mesh_res      = [200 0 0];

%% prepare simulation folder
Sim_Path = 'tmp';
Sim_CSX = 'msl2.xml';
if ~postproc_only
    [status, message, messageid] = rmdir( Sim_Path, 's' ); % clear previous directory
    [status, message, messageid] = mkdir( Sim_Path ); % create empty simulation folder
end

%% setup FDTD parameters & excitation function %%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_timesteps = 20000;
min_decrement = 1e-6;
f_max         = 7e9;
FDTD = InitFDTD( max_timesteps, min_decrement, 'OverSampling', 10 );
FDTD = SetGaussExcite( FDTD, f_max/2, f_max/2 );
BC   = {'MUR' 'MUR' 'PEC' 'PEC' 'PEC' 'PMC'};
FDTD = SetBoundaryCond( FDTD, BC );

%% setup CSXCAD geometry & mesh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CSX = InitCSX();
resolution = c0/(f_max*sqrt(substrate_epr))/unit /50; % resolution of lambda/50
mesh.x = SmoothMeshLines( [-MSL_length MSL_length], resolution );
mesh.y = SmoothMeshLines( [-4*MSL_width -MSL_width/2 MSL_width/2 4*MSL_width], resolution );
mesh.z = SmoothMeshLines( [linspace(0,substrate_thickness,5) 10*substrate_thickness], resolution );
CSX = DefineRectGrid( CSX, unit, mesh );

%% substrate
CSX = AddMaterial( CSX, 'RO4350B' );
CSX = SetMaterialProperty( CSX, 'RO4350B', 'Epsilon', substrate_epr );
start = [mesh.x(1),   mesh.y(1),   0];
stop  = [mesh.x(end), mesh.y(end), substrate_thickness];
CSX = AddBox( CSX, 'RO4350B', 0, start, stop );

%% MSL port
CSX = AddMetal( CSX, 'PEC' );
portstart = [          0, -MSL_width/2, substrate_thickness];
portstop  = [ MSL_length,  MSL_width/2, 0];
[CSX,portstruct] = AddMSLPort( CSX, 999, 1, 'PEC', portstart, portstop, [1 0 0], [0 0 1], [], 'excite' );

%% MSL
start = [-MSL_length, -MSL_width/2, substrate_thickness];
stop  = [          0,  MSL_width/2, substrate_thickness];
CSX = AddBox( CSX, 'PEC', 999, start, stop ); % priority needs to be higher than 

%% define dump boxes
start = [mesh.x(1),   mesh.y(1),   substrate_thickness/2];
stop  = [mesh.x(end), mesh.y(end), substrate_thickness/2];
CSX = AddDump( CSX, 'Et_', 'DumpType', 0,'DumpMode', 2 ); % cell interpolated
CSX = AddBox( CSX, 'Et_', 0, start, stop );
CSX = AddDump( CSX, 'Ht_', 'DumpType', 1,'DumpMode', 2 ); % cell interpolated
CSX = AddBox( CSX, 'Ht_', 0, start, stop );
 
%% write openEMS compatible xml-file
WriteOpenEMS( [Sim_Path '/' Sim_CSX], FDTD, CSX );

%% show the structure
if ~postproc_only
    CSXGeomPlot( [Sim_Path '/' Sim_CSX] );
end

%% run openEMS
openEMS_opts = '';
openEMS_opts = [openEMS_opts ' --engine=fastest'];
% openEMS_opts = [openEMS_opts ' --debug-material'];
% openEMS_opts = [openEMS_opts ' --debug-boxes'];
% openEMS_opts = [openEMS_opts ' --debug-PEC'];
if ~postproc_only
    RunOpenEMS( Sim_Path, Sim_CSX, openEMS_opts );
end


%% postprocess
f = linspace( 1e6, f_max, 1601 );
U = ReadUI( {'port_ut1A','port_ut1B','port_ut1C','et'}, 'tmp/', f );
I = ReadUI( {'port_it1A','port_it1B'}, 'tmp/', f );

% Z = (U.FD{1}.val+U.FD{2}.val)/2 ./ I.FD{1}.val;
% plot( f*1e-9, [real(Z);imag(Z)],'Linewidth',2);
% xlabel('frequency (GHz)');
% ylabel('impedance (Ohm)');
% grid on;
% legend( {'real','imaginary'}, 'location', 'northwest' )
% title( 'line impedance (will fail in case of reflections!)' );

figure
ax = plotyy( U.TD{1}.t/1e-6, [U.TD{1}.val;U.TD{2}.val;U.TD{3}.val], U.TD{4}.t/1e-6, U.TD{4}.val );
xlabel( 'time (us)' );
ylabel( 'amplitude (V)' );
grid on
title( 'Time domain voltage probes and excitation signal' );
legend( {'ut1A','ut1B','ut1C','excitation'} );
% now make the y-axis symmetric to y=0 (align zeros of y1 and y2)
y1 = ylim(ax(1));
y2 = ylim(ax(2));
ylim( ax(1), [-max(abs(y1)) max(abs(y1))] );
ylim( ax(2), [-max(abs(y2)) max(abs(y2))] );

figure
plot( I.TD{1}.t/1e-6, [I.TD{1}.val;I.TD{2}.val] );
xlabel( 'time (us)' );
ylabel( 'amplitude (A)' );
grid on
title( 'Time domain current probes' );
legend( {'it1A','it1B'} );

figure
ax = plotyy( U.FD{1}.f/1e9, abs([U.FD{1}.val;U.FD{2}.val;U.FD{3}.val]), U.FD{1}.f/1e9, angle([U.FD{1}.val;U.FD{2}.val;U.FD{3}.val])/pi*180 );
xlabel( 'frequency (GHz)' );
ylabel( ax(1), 'amplitude (A)' );
ylabel( ax(2), 'phase (deg)' );
grid on
title( 'Frequency domain voltage probes' );
legend( {'abs(uf1A)','abs(uf1B)','abs(uf1C)','angle(uf1A)','angle(uf1B)','angle(uf1C)'} );

figure
ax = plotyy( I.FD{1}.f/1e9, abs([I.FD{1}.val;I.FD{2}.val]), I.FD{1}.f/1e9, angle([I.FD{1}.val;I.FD{2}.val])/pi*180 );
xlabel( 'frequency (GHz)' );
ylabel( ax(1), 'amplitude (A)' );
ylabel( ax(2), 'phase (deg)' );
grid on
title( 'Frequency domain current probes' );
legend( {'abs(if1A)','abs(if1B)','angle(if1A)','angle(if1B)'} );

%% port analysis
[U,I,beta,ZL] = calcPort( portstruct, Sim_Path, f );
%% attention! the reflection coefficient S11 is normalized to ZL!

figure
plot( sin(0:0.01:2*pi), cos(0:0.01:2*pi), 'Color', [.7 .7 .7] );
hold on
plot( 0.5+0.5*sin(0:0.01:2*pi), 0.5*cos(0:0.01:2*pi), 'Color', [.7 .7 .7] );
plot( [-1 1], [0 0], 'Color', [.7 .7 .7] );
plot( S11, 'k' );
plot( real(S11(1)), imag(S11(1)), '*r' );
axis equal
title( 'Reflection coefficient S11 at the measurement plane' );

figure
plot( sin(0:0.01:2*pi), cos(0:0.01:2*pi), 'Color', [.7 .7 .7] );
hold on
plot( 0.5+0.5*sin(0:0.01:2*pi), 0.5*cos(0:0.01:2*pi), 'Color', [.7 .7 .7] );
plot( [-1 1], [0 0], 'Color', [.7 .7 .7] );
Z = vi.FD.v.val ./ vi.FD.i.val;
S11_ = (Z-ZL) ./ (Z+ZL);
plot( S11_, 'k' );
plot( real(S11_(1)), imag(S11_(1)), '*r' );
axis equal
title( {'Reflection coefficient S11 at the measurement plane' 'calculated from voltages and currents'} );

figure
plot( f/1e9, [real(S11);imag(S11)], 'Linewidth',2 );
legend( {'Re(S11)', 'Im(S11)'} );
ylabel( 'amplitude' );
xlabel( 'frequency (GHz)' );
title( 'Reflection coefficient S11 at the measurement plane' );

figure
plotyy( f/1e9, 20*log10(abs(S11)), f/1e9, angle(S11)/pi*180 );
legend( {'|S11|', 'angle(S11)'} );
xlabel( 'frequency (GHz)' );
ylabel( '|S11| (dB)' );
title( 'Reflection coefficient S11 at the measurement plane' );

figure
plot( f/1e9, [real(beta);imag(beta)], 'Linewidth',2 );
legend( 'Re(beta)', 'Im(beta)' );
ylabel( 'propagation constant beta (1/m)' );
xlabel( 'frequency (GHz)' );
title( 'Propagation constant of the MSL' );

figure
plot( f/1e9, [real(ZL);imag(ZL)], 'Linewidth',2);
xlabel('frequency (GHz)');
ylabel('impedance (Ohm)');
grid on;
legend( {'real','imaginary'}, 'location', 'northeast' )
title( 'Characteristic line impedance ZL' );

%% reference plane shift (to the end of the port)
ref_shift = abs(portstop(1) - portstart(1));
[U, I,beta,ZL] = calcPort( portstruct, Sim_Path, f );
%%

figure
plotyy( f/1e9, 20*log10(abs(S11)), f/1e9, angle(S11)/pi*180 );
legend( {'abs(S11)', 'angle(S11)'} );
xlabel( 'frequency (GHz)' );
title( 'Reflection coefficient S11 at the reference plane (at the electric wall)' );

figure
plot( sin(0:0.01:2*pi), cos(0:0.01:2*pi), 'Color', [.7 .7 .7] );
hold on
plot( 0.5+0.5*sin(0:0.01:2*pi), 0.5*cos(0:0.01:2*pi), 'Color', [.7 .7 .7] );
plot( [-1 1], [0 0], 'Color', [.7 .7 .7] );
plot( S11, 'k' );
plot( real(S11(1)), imag(S11(1)), '*r' );
axis equal
title( 'Reflection coefficient S11 at the reference plane (at the electric wall)' );

figure
plot( sin(0:0.01:2*pi), cos(0:0.01:2*pi), 'Color', [.7 .7 .7] );
hold on
plot( 0.5+0.5*sin(0:0.01:2*pi), 0.5*cos(0:0.01:2*pi), 'Color', [.7 .7 .7] );
plot( [-1 1], [0 0], 'Color', [.7 .7 .7] );
Z = vi.FD.v.val_shifted ./ vi.FD.i.val_shifted;
S11_ = (Z-ZL) ./ (Z+ZL);
plot( S11_, 'k' );
plot( real(S11_(1)), imag(S11_(1)), '*r' );
axis equal
title( {'Reflection coefficient S11 at the reference plane (at the electric wall)' 'calculated from shifted voltages and currents'} );

%% visualize electric and magnetic fields
% you will find vtk dump files in the simulation folder (tmp/)
% use paraview to visualize them
