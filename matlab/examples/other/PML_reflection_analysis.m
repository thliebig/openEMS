%
% fake-PML parallel plate waveguide example
%
% this example analyzes the reflection coefficient of a vacuum-pml
% interface
%

%
% currently this example uses a normal material with a certain conductivity
% profile and not a pml
%

close all
% clear
clc

physical_constants


postprocessing_only = 0;



%% setup the simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawingunit = 1e-6; % specify everything in um

length = 10000;
epr = 1;

mesh_res = [200 200 200];
max_timesteps = 100000;
min_decrement = 1e-6;
f_max = 8e9;

%% setup FDTD parameters & excitation function %%%%%%%%%%%%%%%%%%%%%%%%%%%%
FDTD = InitFDTD( max_timesteps, min_decrement );
FDTD = SetGaussExcite( FDTD, f_max/2, f_max/2 );
BC = [0 0 1 1 0 0];
FDTD = SetBoundaryCond( FDTD, BC );

%% mesh grading
N_pml = 8;
pml_delta = cumsum(mesh_res(1) * 1.0 .^ (1:N_pml));
% pml_delta = cumsum([200 200 200 200 200]);

%% setup CSXCAD geometry & mesh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CSX = InitCSX();
mesh.x = 0 : mesh_res(1) : length;
mesh.x = [mesh.x(1) - fliplr(pml_delta), mesh.x];
mesh.y = -2*mesh_res(2) : mesh_res(2) : 2*mesh_res(2);
mesh.z = 0 : mesh_res(3) : 4*mesh_res(3);
CSX = DefineRectGrid( CSX, drawingunit, mesh );

%% fake pml %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g = 2; % 2..3
R0 = 1e-6; % requested analytical reflection coefficient
Zm = sqrt(MUE0/(EPS0*epr)); % calculate reflection for substrate/pml interface
delta = pml_delta(end) * drawingunit;
deltal = mean(diff(pml_delta)) * drawingunit;
kappa0 = -log(R0)*log(g)/( 2*Zm*deltal*(g^(delta/deltal)-1) );

% kappa0 = 1.05;
CSX = AddMaterial( CSX, 'pml_xmin' );
CSX = SetMaterialProperty( CSX, 'pml_xmin', 'Epsilon', epr );
CSX = SetMaterialProperty( CSX, 'pml_xmin', 'Kappa', kappa0 );
CSX = SetMaterialProperty( CSX, 'pml_xmin', 'Sigma', kappa0 * MUE0/(EPS0*epr) );
CSX = SetMaterialWeight( CSX, 'pml_xmin', 'Kappa', [num2str(g) '^((abs(x-100)-' num2str(abs(mesh.x(N_pml+1))) ')/(' num2str(deltal) '/' num2str(drawingunit) '))'] ); % g^(rho/deltal)*kappa0
CSX = SetMaterialWeight( CSX, 'pml_xmin', 'Sigma', [num2str(g) '^((abs(x-100)-' num2str(abs(mesh.x(N_pml+1))) ')/(' num2str(deltal) '/' num2str(drawingunit) '))'] );
start = [mesh.x(1), mesh.y(1),   mesh.z(1)];
stop  = [100, mesh.y(end), mesh.z(end)];
CSX = AddBox( CSX, 'pml_xmin', 1, start, stop );

figure
x = [-fliplr(pml_delta) 50];
plot( x, kappa0 * g.^((abs(x-50)-abs(mesh.x(N_pml+1)))./(deltal/drawingunit)) ,'x-');
xlabel( 'x / m' );
ylabel( 'kappa' );
figure
title( 'conductivity profile inside the material' );

%% excitation
CSX = AddExcitation( CSX, 'excitation1', 0, [0 0 1]);
idx = interp1( mesh.x, 1:numel(mesh.x), length*2/3, 'nearest' );
start = [mesh.x(idx), mesh.y(1),   mesh.z(1)];
stop  = [mesh.x(idx), mesh.y(end), mesh.z(end)];
CSX = AddBox( CSX, 'excitation1', 0, start, stop );

%% define dump boxes... %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CSX = AddDump( CSX, 'Et_', 'DumpMode', 2 );
start = [mesh.x(1),   mesh.y(1),   mesh.z(3)];
stop  = [mesh.x(end), mesh.y(end), mesh.z(3)];
CSX = AddBox( CSX, 'Et_', 0, start, stop );

CSX = AddDump( CSX, 'Ht_', 'DumpType', 1, 'DumpMode', 2 );
CSX = AddBox( CSX, 'Ht_', 0, start, stop );

% hdf5 file
CSX = AddDump( CSX, 'E', 'DumpType', 0, 'DumpMode', 2, 'FileType', 1 );
idx = interp1( mesh.x, 1:numel(mesh.x), length*1/3, 'nearest' );
start = [mesh.x(idx), mesh.y(3), mesh.z(1)];
stop  = [mesh.x(idx), mesh.y(3), mesh.z(end)];
CSX = AddBox( CSX, 'E', 0, start, stop );

% hdf5 file
CSX = AddDump( CSX, 'H', 'DumpType', 1, 'DumpMode', 2, 'FileType', 1 );
idx = interp1( mesh.x, 1:numel(mesh.x), length*1/3, 'nearest' );
start = [mesh.x(idx), mesh.y(1),   mesh.z(3)];
stop  = [mesh.x(idx), mesh.y(end), mesh.z(3)];
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
Sim_CSX = 'PML_reflection_analysis.xml';

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
f = linspace( 0, f_max, 1601 );
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


Z0 = sqrt(MUE0/EPS0); % line impedance
Z = Ef ./ Hf; % impedance at measurement plane
gamma = (Z - Z0) ./ (Z + Z0);

plot( f/1e9, 20*log10(abs(gamma)),'Linewidth',2);
xlabel('frequency (GHz)');
ylabel('reflection coefficient gamma (dB)');
grid on;
title( 'Reflection Coefficient' );

if exist('ref_1','var')
    hold on
    plot( f/1e9, ref_1,'--','Linewidth',2, 'Color', [1 0 0]);
    hold off
end
ref_1 = 20*log10(abs(gamma));
