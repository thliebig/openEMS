%
% Tutorials / 7T MRI Loop Coil
%
% Description at:
% http://openems.de/index.php/Tutorial:_MRI_Loop_Coil
%
% Tested with
%  - openEMS v0.0.33
%  - Matlab 7.12.0 (R2011a)
%
% (C) 2013-2015 Thorsten Liebig <thorsten.liebig@gmx.de>

close all
clear
clc

%% setup the simulation
physical_constants; %get some physical constants like c0 and MUE0
unit = 1e-3; % all length in mm

% Loop-Coil parameter
loop.length = 80;       % length of the loop (in z-direction)
loop.width = 60;        % width of the loop  (in y-direction)
loop.strip_width = 5;   % metal strip width
loop.strip_N_cells = 3; % number of cells over the strip length
loop.air_gap = loop.strip_width/3;       % air gap width for lumped capacitors
loop.pos_x = -130;       % position of loop
loop.C_gap = 5.4e-12;   % lumped cap value
loop.port_R = 10;       % feeding port resistance

%% define the human body model (virtual family)
% set file name for human body model to create with "Convert_VF_DiscMaterial"
% the file name should contain a full path
body_model_file = [pwd '/Ella_centered_298MHz.h5'];

% convert only part of the model (head/shoulder section)
body_model_range = {[],[],[-0.85 -0.4]};

% paths to virtual family voxel models (VFVM), adept to your install!
VF_raw_filesuffix = '/tmp/Ella_26y_V2_1mm';
VF_mat_db_file = '/tmp/DB_h5_20120711_SEMCADv14.8.h5';

% delete(body_model_file); % uncomment to delete old model if something changed

% convert model (if it does not exist)
Convert_VF_DiscMaterial(VF_raw_filesuffix, VF_mat_db_file, body_model_file, ...
                        'Frequency', 298e6, 'Center', 1, ...
                        'Range', body_model_range);

% rotate model to face the nose in x-dir, and translate
body_model_transform = {'Rotate_X',pi,'Rotate_Z',pi/2, ...
                        'Translate',[0,5,-720]};

% the head should + part of shoulder should fit this box
body_box.start = [-120 -150 -200];
body_box.stop  = [+100 +150 +130];

% box with high res mesh
mesh_box.start = [-120 -80 -120];
mesh_box.stop  = [+100 +80 +120];
mesh_box.resolution = 2;

%% some mesh parameter
Air_Box = 150;      % size of the surrounding air box (150mm)

%% setup FDTD parameter & excitation function
% init FDTD structure
FDTD = InitFDTD( 'EndCriteria', 1e-4, 'CellConstantMaterial', 0);

% define gaussian pulse excitation signal
f0 = 298e6; % center frequency
fc = 300e6; % 20 dB corner frequency
FDTD = SetGaussExcite( FDTD, f0, fc );

% setup boundary conditions
BC = {'MUR' 'MUR' 'MUR' 'MUR' 'MUR' 'MUR'}; % boundary conditions
FDTD = SetBoundaryCond( FDTD, BC );

%% setup CSXCAD geometry & mesh
CSX = InitCSX();

%% create loop
% setup all properties needed
CSX = AddMetal( CSX, 'loop' );
CSX = AddLumpedElement( CSX, 'caps_y', 1, 'C', loop.C_gap);
CSX = AddLumpedElement( CSX, 'caps_z', 2, 'C', loop.C_gap);

% horizontal (y-direction) strips
start = [loop.pos_x -loop.width/2   -loop.length/2];
stop  = [loop.pos_x -loop.air_gap/2 -loop.length/2+loop.strip_width];
CSX = AddBox(CSX,'loop',10,start,stop);

start = [loop.pos_x  -loop.width/2   loop.length/2 ];
stop  = [loop.pos_x  -loop.air_gap/2 loop.length/2-loop.strip_width];
CSX = AddBox(CSX,'loop',10,start,stop);

start = [loop.pos_x loop.width/2   -loop.length/2];
stop  = [loop.pos_x loop.air_gap/2 -loop.length/2+loop.strip_width];
CSX = AddBox(CSX,'loop',10,start,stop);

start = [loop.pos_x loop.width/2   loop.length/2  ];
stop  = [loop.pos_x loop.air_gap/2 loop.length/2-loop.strip_width];
CSX = AddBox(CSX,'loop',10,start,stop);

% vertical (z-direction) strips
start = [loop.pos_x -loop.width/2                  -loop.length/2+loop.strip_width];
stop  = [loop.pos_x -loop.width/2+loop.strip_width -loop.air_gap/2];
CSX = AddBox(CSX,'loop',10,start,stop);

start = [loop.pos_x -loop.width/2                  loop.length/2-loop.strip_width];
stop  = [loop.pos_x -loop.width/2+loop.strip_width loop.air_gap/2];
CSX = AddBox(CSX,'loop',10,start,stop);

start = [loop.pos_x loop.width/2                   -loop.length/2+loop.strip_width];
stop  = [loop.pos_x loop.width/2-loop.strip_width  -loop.air_gap/2];
CSX = AddBox(CSX,'loop',10,start,stop);

start = [loop.pos_x loop.width/2                   loop.length/2-loop.strip_width ];
stop  = [loop.pos_x loop.width/2-loop.strip_width  loop.air_gap/2];
CSX = AddBox(CSX,'loop',10,start,stop);

% add the lumped capacities
start = [loop.pos_x -loop.width/2+loop.strip_width/2-loop.air_gap/2 -loop.air_gap/2];
stop  = [loop.pos_x -loop.width/2+loop.strip_width/2+loop.air_gap/2 +loop.air_gap/2];
CSX = AddBox(CSX,'caps_z',10,start,stop);

start = [loop.pos_x loop.width/2-loop.strip_width/2-loop.air_gap/2 -loop.air_gap/2];
stop  = [loop.pos_x loop.width/2-loop.strip_width/2+loop.air_gap/2 +loop.air_gap/2];
CSX = AddBox(CSX,'caps_z',10,start,stop);

start = [loop.pos_x -loop.air_gap/2 loop.length/2-loop.strip_width/2-loop.air_gap/2];
stop  = [loop.pos_x +loop.air_gap/2 loop.length/2-loop.strip_width/2+loop.air_gap/2];
CSX = AddBox(CSX,'caps_y',10,start,stop);

% add a lumped port as excitation
start = [loop.pos_x -loop.air_gap/2 -loop.length/2+loop.strip_width/2-loop.air_gap/2];
stop  = [loop.pos_x +loop.air_gap/2 -loop.length/2+loop.strip_width/2+loop.air_gap/2];
[CSX port] = AddLumpedPort(CSX, 100, 1, loop.port_R, start, stop, [0 1 0], true);

%% define human body model
CSX = AddDiscMaterial(CSX, 'body_model', 'File', body_model_file, 'Scale', 1/unit, 'Transform', body_model_transform);
CSX = AddBox(CSX, 'body_model', 0, body_box.start, body_box.stop);

%% finalize mesh
% create loop mesh
mesh = DetectEdges(CSX);

% add a dense homogeneous mesh inside the human body model
mesh.x = [mesh.x mesh_box.start(1) mesh_box.stop(1)];
mesh.y = [mesh.y mesh_box.start(2) mesh_box.stop(2)];
mesh.z = [mesh.z mesh_box.start(3) mesh_box.stop(3)];

% add lines in x-dir for the loop and a cell centered around 0
mesh.x = [mesh.x loop.pos_x -mesh_box.resolution/2 mesh_box.resolution/2];

% smooth the mesh for the loop & body
mesh = SmoothMesh(mesh, mesh_box.resolution);

% add air spacer
mesh.x = [-Air_Box+mesh.x(1) mesh.x mesh.x(end)+Air_Box];
mesh.y = [-Air_Box+mesh.y(1) mesh.y mesh.y(end)+Air_Box];
mesh.z = [-Air_Box+mesh.z(1) mesh.z mesh.z(end)+Air_Box];

mesh = SmoothMesh(mesh, c0 / (f0+fc) / unit / 10, 1.5, 'algorithm', 1);

%% Add Dump boxes (2D boxes) for H and SAR on xy- and xz-plane
CSX = AddDump(CSX,'Hf_xy','DumpType',11,'FileType',1,'Frequency',f0);
CSX = AddBox(CSX,'Hf_xy',0, body_box.start.*[1 1 0], body_box.stop.*[1 1 0]);
CSX = AddDump(CSX,'SAR_xy','DumpType',20,'DumpMode',2,'FileType',1,'Frequency',f0);
CSX = AddBox(CSX,'SAR_xy',0, body_box.start.*[1 1 0], body_box.stop.*[1 1 0]);

CSX = AddDump(CSX,'Hf_xz','DumpType',11,'FileType',1,'Frequency',f0);
CSX = AddBox(CSX,'Hf_xz',0, body_box.start.*[1 0 1], body_box.stop.*[1 0 1]);
CSX = AddDump(CSX,'SAR_xz','DumpType',20,'DumpMode',2,'FileType',1,'Frequency',f0);
CSX = AddBox(CSX,'SAR_xz',0, body_box.start.*[1 0 1], body_box.stop.*[1 0 1]);

%% add 10 lines in all direction to make space for PML or MUR absorbing
%% boundary conditions
mesh = AddPML(mesh, 10);

%% finally define the FDTD mesh grid
disp(['number of cells: ' num2str(1e-6*numel(mesh.x)*numel(mesh.y)*numel(mesh.z)) ' Mcells'])
CSX = DefineRectGrid( CSX, unit, mesh );

%% prepare simulation folder
Sim_Path = ['tmp_' mfilename];
Sim_CSX = [mfilename '.xml'];

[status, message, messageid] = rmdir( Sim_Path, 's' ); % clear previous directory
[status, message, messageid] = mkdir( Sim_Path ); % create empty simulation folder

%% write openEMS compatible xml-file
WriteOpenEMS( [Sim_Path '/' Sim_CSX], FDTD, CSX );

%% show the structure and export as vtk data automatically
CSXGeomPlot( [Sim_Path '/' Sim_CSX] , ['--export-polydata-vtk=' Sim_Path ' --RenderDiscMaterial -v']);

%% run openEMS
RunOpenEMS( Sim_Path, Sim_CSX);

%% postprocessing & do the plots
freq = linspace( f0-fc, f0+fc, 501 );
port = calcPort(port, Sim_Path, freq);

Zin = port.uf.tot ./ port.if.tot;
s11 = port.uf.ref ./ port.uf.inc;

% get the feeding power for frequency f0
P0_in = interp1(freq, port.P_acc, f0);

%%
% plot reflection coefficient S11
figure
h = plot( freq/1e6, 20*log10(abs(s11)), 'k-', 'Linewidth', 2 );
grid on
title( 'reflection coefficient S_{11}' );
xlabel( 'frequency f / MHz' );
ylabel( 'reflection coefficient |S_{11}| (dB)' );

% plot feed point admittance
figure
h = plot( freq/1e6, real(1./Zin), 'k-', 'Linewidth', 2 );
hold on
grid on
plot( freq/1e6, imag(1./Zin), 'r--', 'Linewidth', 2 );
title( 'feed port admittance' );
xlabel( 'frequency f (MHz)' );
ylabel( 'admittance Y_{in} (S)' );
legend( 'real', 'imag' );

%% read SAR values on a xy-plane (range)
[SAR SAR_mesh] = ReadHDF5Dump([Sim_Path '/SAR_xy.h5']);
SAR = SAR.FD.values{1}/P0_in;

% SAR plot
figure()
subplot(1,2,1);
[X Y] = ndgrid(SAR_mesh.lines{1},SAR_mesh.lines{2});
colormap('hot');
h = pcolor(X,Y,(squeeze(SAR)));
% h = pcolor(X,Y,log10(squeeze(SAR)));
set(h,'EdgeColor','none');
xlabel('x -->');
ylabel('y -->');
title('local SAR');
axis equal tight

%% read SAR values on a xz-plane (range)
[SAR SAR_mesh] = ReadHDF5Dump([Sim_Path '/SAR_xz.h5']);
SAR = SAR.FD.values{1}/P0_in;

% SAR plot
subplot(1,2,2);
[X Z] = ndgrid(SAR_mesh.lines{1},SAR_mesh.lines{3});
colormap('hot');
h = pcolor(X,Z,(squeeze(SAR)));
% h = pcolor(X,Y,log10(squeeze(SAR)));
set(h,'EdgeColor','none');
xlabel('x -->');
ylabel('z -->');
title('local SAR');
axis equal tight

%% plot B1+/- on an xy-plane
[H_field H_mesh] = ReadHDF5Dump([Sim_Path '/Hf_xy.h5']);
% calc Bx,By, B1p, B1m normalize to the input-power
Bx = MUE0*H_field.FD.values{1}(:,:,:,1)/sqrt(P0_in);
By = MUE0*H_field.FD.values{1}(:,:,:,2)/sqrt(P0_in);
B1p = 0.5*(Bx+1j*By);
B1m = 0.5*(Bx-1j*By);
% create a 2D grid to plot on
[X Y] = ndgrid(H_mesh.lines{1},H_mesh.lines{2});

Dump2VTK([Sim_Path '/B1p_xy.vtk'], abs(B1p), H_mesh, 'B-Field');
Dump2VTK([Sim_Path '/B1m_xy.vtk'], abs(B1m), H_mesh, 'B-Field');

% B1+ plot
figure()
subplot(1,2,1);
h = pcolor(X,Y,log10(abs(B1p)));
set(h,'EdgeColor','none');
xlabel('x -->');
ylabel('y -->');
title('B_1^+ field (dB)');
axis equal tight

% B1- plot
subplot(1,2,2);
h = pcolor(X,Y,log10(abs(B1m)));
set(h,'EdgeColor','none');
xlabel('x -->');
ylabel('y -->');
title('B_1^- field (dB)');
axis equal tight

%% plot B1+/- on an xz-plane
[H_field H_mesh] = ReadHDF5Dump([Sim_Path '/Hf_xz.h5']);
% calc Bx,By, B1p, B1m normalize to the input-power
Bx = MUE0*H_field.FD.values{1}(:,:,:,1)/sqrt(P0_in);
By = MUE0*H_field.FD.values{1}(:,:,:,2)/sqrt(P0_in);
B1p = 0.5*(Bx+1j*By);
B1m = 0.5*(Bx-1j*By);
% create a 2D grid to plot on
[X Z] = ndgrid(H_mesh.lines{1},H_mesh.lines{3});

Dump2VTK([Sim_Path '/B1p_xz.vtk'], abs(B1p), H_mesh, 'B-Field');
Dump2VTK([Sim_Path '/B1m_xz.vtk'], abs(B1m), H_mesh, 'B-Field');

% B1+ plot
figure()
subplot(1,2,1);
h = pcolor(X,Z,log10(squeeze(abs(B1p))));
set(h,'EdgeColor','none');
xlabel('x -->');
ylabel('z -->');
title('B_1^+ field (dB)');
axis equal tight

% B1- plot
subplot(1,2,2);
h = pcolor(X,Z,log10(squeeze(abs(B1m))));
set(h,'EdgeColor','none');
xlabel('x -->');
ylabel('z -->');
title('B_1^- field (dB)');
axis equal tight

%% dump to vtk to view in Paraview
ConvertHDF5_VTK([Sim_Path '/SAR_xy.h5'],[Sim_Path '/SAR_xy'], 'weight', 1/P0_in, 'FieldName', 'SAR');
ConvertHDF5_VTK([Sim_Path '/SAR_xz.h5'],[Sim_Path '/SAR_xz'], 'weight', 1/P0_in, 'FieldName', 'SAR');

%%
ConvertHDF5_VTK([Sim_Path '/Hf_xy.h5'],[Sim_Path '/B1_xy'], 'weight', MUE0/sqrt(P0_in), 'FieldName', 'B1-field');
ConvertHDF5_VTK([Sim_Path '/Hf_xz.h5'],[Sim_Path '/B1_xz'], 'weight', MUE0/sqrt(P0_in), 'FieldName', 'B1-field');
