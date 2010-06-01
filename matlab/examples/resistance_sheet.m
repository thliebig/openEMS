%
% resistance "sheet" example
%
% this example calculates the reflection coefficient of a sheet resistance
% at the end of a parallel plate wave guide
%
% play around with the R and epr values
%

close all
clear
clc

physical_constants


postprocessing_only = 0;



%% setup the simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
epr = 1; % relative permittivity of the material inside the parallel plate waveguide

% define the resistance
R = sqrt(MUE0/(EPS0*epr)); % matched load (no reflections) (vacuum: approx. 377 Ohm)
% R = 1e-10; % short circuit (reflection coefficient = -1)
% R = 1e10; % open circuit (reflection coefficient = 1)


drawingunit = 1e-6; % specify everything in um
length = 10000;
mesh_res = [200 200 200];
max_timesteps = 100000;
min_decrement = 1e-6;
f_max = 1e9;

%% setup FDTD parameters & excitation function %%%%%%%%%%%%%%%%%%%%%%%%%%%%
FDTD = InitFDTD( max_timesteps, min_decrement );
FDTD = SetGaussExcite( FDTD, f_max/2, f_max/2 );
BC = [1 2 1 1 0 0]; % 0:PEC  1:PMC  2:MUR-ABC
FDTD = SetBoundaryCond( FDTD, BC );

%% setup CSXCAD geometry & mesh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CSX = InitCSX();
mesh.x = 0 : mesh_res(1) : length;
mesh.y = -2*mesh_res(2) : mesh_res(2) : 2*mesh_res(2);
mesh.z = 0 : mesh_res(3) : 4*mesh_res(3);
CSX = DefineRectGrid( CSX, drawingunit, mesh );

%% measurement plane & reference plane
meas_plane_xidx = interp1( mesh.x, 1:numel(mesh.x), length*1/3, 'nearest' );
ref_plane_xidx = 3;

%% fill the parallel plate waveguide with material
CSX = AddMaterial( CSX, 'm1' );
CSX = SetMaterialProperty( CSX, 'm1', 'Epsilon', epr );
start = [mesh.x(1),   mesh.y(1),   mesh.z(1)];
stop  = [mesh.x(end), mesh.y(end), mesh.z(end)];
CSX = AddBox( CSX, 'm1', -1, start, stop );

%% excitation
CSX = AddExcitation( CSX, 'excitation1', 0, [0 0 1]);
idx = interp1( mesh.x, 1:numel(mesh.x), length*2/3, 'nearest' );
start = [mesh.x(idx), mesh.y(1),   mesh.z(1)];
stop  = [mesh.x(idx), mesh.y(end), mesh.z(end)];
CSX = AddBox( CSX, 'excitation1', 0, start, stop );

%% define the sheet resistance
start = [mesh.x(ref_plane_xidx-1), mesh.y(1),   mesh.z(1)];
stop  = [mesh.x(ref_plane_xidx),   mesh.y(end), mesh.z(end)];
l = abs(mesh.z(end) - mesh.z(1)) * drawingunit; % length of the "sheet"
A = abs(start(1) - stop(1)) * abs(mesh.y(end) - mesh.y(1)) * drawingunit^2; % area of the "sheet"
kappa = l/A / R; % [kappa] = S/m
CSX = AddMaterial( CSX, 'sheet_resistance' );
CSX = SetMaterialProperty( CSX, 'sheet_resistance', 'Kappa', kappa );
CSX = AddBox( CSX, 'sheet_resistance', 0, start, stop );

%% define dump boxes... %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CSX = AddDump( CSX, 'Et_', 'DumpMode', 2 );
start = [mesh.x(1),   mesh.y(1),   mesh.z(3)];
stop  = [mesh.x(end), mesh.y(end), mesh.z(3)];
CSX = AddBox( CSX, 'Et_', 0, start, stop );

CSX = AddDump( CSX, 'Ht_', 'DumpType', 1, 'DumpMode', 2 );
CSX = AddBox( CSX, 'Ht_', 0, start, stop );

% hdf5 file
CSX = AddDump( CSX, 'E', 'DumpType', 0, 'DumpMode', 2, 'FileType', 1 );
start = [mesh.x(meas_plane_xidx), mesh.y(3), mesh.z(1)];
stop  = [mesh.x(meas_plane_xidx), mesh.y(3), mesh.z(end)];
CSX = AddBox( CSX, 'E', 0, start, stop );

% hdf5 file
CSX = AddDump( CSX, 'H', 'DumpType', 1, 'DumpMode', 2, 'FileType', 1 );
start = [mesh.x(meas_plane_xidx), mesh.y(1),   mesh.z(3)];
stop  = [mesh.x(meas_plane_xidx), mesh.y(end), mesh.z(3)];
CSX = AddBox( CSX, 'H', 0, start, stop );

%% define openEMS options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
openEMS_opts = '';
% openEMS_opts = [openEMS_opts ' --disable-dumps'];
% openEMS_opts = [openEMS_opts ' --debug-material'];
% openEMS_opts = [openEMS_opts ' --debug-operator'];
% openEMS_opts = [openEMS_opts ' --debug-boxes'];
% openEMS_opts = [openEMS_opts ' --showProbeDiscretization'];
openEMS_opts = [openEMS_opts ' --engine=fastest'];

Sim_Path = 'tmp';
Sim_CSX = 'tmp.xml';

if ~postprocessing_only
    [~,~,~] = rmdir(Sim_Path,'s');
    [~,~,~] = mkdir(Sim_Path);
end

%% Write openEMS compatible xml-file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WriteOpenEMS([Sim_Path '/' Sim_CSX],FDTD,CSX);

%% cd to working dir and run openEMS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~postprocessing_only
    savePath = pwd;
    cd(Sim_Path); %cd to working dir
    args = [Sim_CSX ' ' openEMS_opts];
    invoke_openEMS(args);
    cd(savePath)
end


%% postproc & do the plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% E_coords = ReadHDF5Mesh( [Sim_Path '/E.h5'] );
% H_coords = ReadHDF5Mesh( [Sim_Path '/H.h5'] );
E = ReadHDF5FieldData( [Sim_Path '/E.h5'] );
H = ReadHDF5FieldData( [Sim_Path '/H.h5'] );
E_val = cellfun( @(x) squeeze(x(1,1,:,3)), E.values, 'UniformOutput', false );
H_val = cellfun( @(x) squeeze(x(1,:,1,2)), H.values, 'UniformOutput', false );
E_val = cell2mat(E_val);
H_val = cell2mat(H_val.');

% pick center point
Et = E_val(3,:);
Ht = H_val(:,3).';

delta_t_2 = H.time(1) - E.time(1); % half time-step (s) 

% create finer frequency resolution
f = linspace( 0, f_max, 201 );
Ef = DFT_time2freq( E.time, Et, f );
Hf = DFT_time2freq( H.time, Ht, f );
Hf = Hf .* exp(-1i*2*pi*f*delta_t_2); % compensate half time-step advance of H-field

% H is now time interpolated, but the position is not corrected with
% respect to E

% figure
% plot( E.time/1e-6, Et );
% xlabel('time (us)');
% ylabel('amplitude (V)');
% grid on;
% title( 'Time domain voltage probe' );
% 
% figure
% plot( H.time/1e-6, Ht );
% xlabel('time (us)');
% ylabel('amplitude (A)');
% grid on;
% title( 'Time domain current probe' );


Z0 = sqrt(MUE0/(EPS0*epr)); % line impedance
Z = Ef ./ Hf; % impedance at measurement plane
gamma = (Z - Z0) ./ (Z + Z0);

% reference plane shift
beta = 2*pi*f * sqrt(MUE0*(EPS0*epr)); % TEM wave
meas_plane_x = mesh.x(meas_plane_xidx);
ref_plane_x = mesh.x(ref_plane_xidx);
gamma_refplane = gamma .* exp(2i*beta* (meas_plane_x-ref_plane_x)*drawingunit);
Z_refplane = Z0 * (1+gamma_refplane)./(1-gamma_refplane);

% smith chart
figure
if exist( 'smith', 'file' )
    % smith chart
    % www.ece.rutgers.edu/~orfanidi/ewa
    % or cmt toolbox from git.ate.uni-duisburg.de
    smith
else
    % poor man smith chart
    plot( sin(0:0.01:2*pi), cos(0:0.01:2*pi), 'Color', [.7 .7 .7] );
    hold on
%     plot( 0.25+0.75*sin(0:0.01:2*pi), 0.75*cos(0:0.01:2*pi), 'Color', [.7 .7 .7] );
    plot( 0.5+0.5*sin(0:0.01:2*pi), 0.5*cos(0:0.01:2*pi), 'Color', [.7 .7 .7] );
%     plot( 0.75+0.25*sin(0:0.01:2*pi), 0.25*cos(0:0.01:2*pi), 'Color', [.7 .7 .7] );
    plot( [-1 1], [0 0], 'Color', [.7 .7 .7] );
    axis equal
end
plot( real(gamma_refplane), imag(gamma_refplane), 'r*' );
% plot( real(gamma), imag(gamma), 'k*' );
title( 'reflection coefficient S11 at reference plane' )

figure
plot( f/1e9, [real(Z_refplane);imag(Z_refplane)],'Linewidth',2);
xlabel('frequency (GHz)');
ylabel('impedance (Ohm)');
grid on;
title( 'Impedance at reference plane' );
legend( {'real','imag'} );
