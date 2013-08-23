%
% Tutorials / 7T MRI Loop Coil
%
% Describtion at:
% http://openems.de/index.php/Tutorial:_MRI_Loop_Coil
%
% Tested with
%  - openEMS v0.0.31
%
% (C) 2013 Thorsten Liebig <thorsten.liebig@uni-due.de>

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
loop.pos_x = 110;       % position of loop
loop.C_gap = 5e-12;     % lumped cap value
loop.port_R = 10;       % feeding port resistance

%% define the phantom material stackup material
materials.name{1}='skin';
materials.rad(1)=100;
materials.eps_r(1)=49.8;
materials.kappa(1)=0.64;
materials.density(1)=1100;
materials.prio(1)=10;

materials.name{2}='headbone';
materials.rad(2)=95;
materials.eps_r(2)=13.4;
materials.kappa(2)=0.08;
materials.density(2)=1990;
materials.prio(2)=11;

materials.name{3}='CSF';
materials.rad(3)=90;
materials.eps_r(3)=72.734;
materials.kappa(3)=2.2245;
materials.density(3)=1007;
materials.prio(3)=12;

materials.name{4}='brain';    % average of white/grey matter
materials.rad(4)=86;
materials.eps_r(4)=60;
materials.kappa(4)=0.69;
materials.density(4)=1039;
materials.prio(4)=13;

%% some mesh parameter
mat_mesh = 2;       % mesh inside the phantom (2mm)
Air_Box = 150;      % size of the surrounding air box (150mm)

%% setup FDTD parameter & excitation function
% init FDTD structure
FDTD = InitFDTD( 'EndCriteria', 1e-4 );

% define gaussian pulse excitation signal
f0 = 300e6; % center frequency
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

%% define materials in a loop
for n=1:numel(materials.name)
    CSX = AddMaterial( CSX, materials.name{n} );
    CSX = SetMaterialProperty( CSX, materials.name{n}, 'Epsilon', materials.eps_r(n), 'Kappa', materials.kappa(n), 'Density', materials.density(n));
    CSX = AddSphere( CSX, materials.name{n}, materials.prio(n), [0 0 0], materials.rad(n),'Transform',{'Scale',[1 0.8 1] } );
end

%% finalize mesh
% create loop mesh (detect only metal and the lumped elements)
mesh = DetectEdges(CSX, [], 'SetPropertyType', {'Metal','LumpedElement'});

% add a dense homegeneous mesh inside the phantom
mesh.x = [mesh.x -materials.rad -mat_mesh/2 mat_mesh/2 materials.rad];
mesh.y = [mesh.y -materials.rad*0.8  materials.rad*0.8];
mesh.z = [mesh.z -materials.rad  materials.rad];

% smooth the mesh for the loop & phantom
mesh = SmoothMesh(mesh, mat_mesh);

% add air spacer
mesh.x = [-Air_Box+mesh.x(1) mesh.x mesh.x(end)+Air_Box];
mesh.y = [-Air_Box+mesh.y(1) mesh.y mesh.y(end)+Air_Box];
mesh.z = [-Air_Box+mesh.z(1) mesh.z mesh.z(end)+Air_Box];

mesh = SmoothMesh(mesh, c0 / (f0+fc) / unit / 10, 1.5, 'algorithm', 1);

%% Add Dump boxes (3D box) for E,J and SAR
CSX = AddDump(CSX,'Hf','DumpType',11,'FileType',1,'Frequency',f0);
CSX = AddDump(CSX,'SAR','DumpType',20,'DumpMode',2,'FileType',1,'Frequency',f0);
start = [-120 -120 -120];
stop  = [+120 +120 +120];
CSX = AddBox(CSX,'Hf',0,start,stop);
CSX = AddBox(CSX,'SAR',0,start,stop);

%% add a nf2ff calc box; size is 3 cells away from MUR boundary condition
start = [mesh.x(1)   mesh.y(1)   mesh.z(1)];
stop  = [mesh.x(end) mesh.y(end) mesh.z(end)];
[CSX nf2ff] = CreateNF2FFBox(CSX, 'nf2ff', start, stop,'Frequency',f0,'OptResolution',c0/f0/unit/20);

%% add 10 lines in all direction to make space for PML or MUR absorbing
%% boundary conditions
mesh = AddPML(mesh, 10);

%% finaly define the FDTD mesh grid
disp(['number of cells: ' num2str(1e-6*numel(mesh.x)*numel(mesh.y)*numel(mesh.z)) ' Mcells'])
CSX = DefineRectGrid( CSX, unit, mesh );

%% prepare simulation folder
Sim_Path = mfilename;
Sim_CSX = [Sim_Path '.xml'];

[status, message, messageid] = rmdir( Sim_Path, 's' ); % clear previous directory
[status, message, messageid] = mkdir( Sim_Path ); % create empty simulation folder

%% write openEMS compatible xml-file
WriteOpenEMS( [Sim_Path '/' Sim_CSX], FDTD, CSX );

%% show the structure and export as vtk data automatically
CSXGeomPlot( [Sim_Path '/' Sim_CSX] , ['--export-polydata-vtk=' Sim_Path]);

%% run openEMS
RunOpenEMS( Sim_Path, Sim_CSX);

%% postprocessing & do the plots
freq = linspace( f0-fc, f0+fc, 501 );
port = calcPort(port, Sim_Path, freq);

Zin = port.uf.tot ./ port.if.tot;
s11 = port.uf.ref ./ port.uf.inc;

P_in = real(0.5 * port.uf.tot .* conj(port.if.tot)); % antenna feed power
% get the feeding power for frequency f0
P0_in = interp1(freq, P_in, f0);

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
[SAR SAR_mesh] = ReadHDF5Dump([Sim_Path '/SAR.h5'], 'Range', {[],[],0});
SAR = SAR.FD.values{1};

%%
% SAR plot
figure()
[X Y] = ndgrid(SAR_mesh.lines{1},SAR_mesh.lines{2});
colormap('hot');
h = pcolor(X,Y,(squeeze(SAR)));
% h = pcolor(X,Y,log10(squeeze(SAR)));
set(h,'EdgeColor','none');
xlabel('x -->');
ylabel('y -->');
title('local SAR');
axis equal tight

%%
[H_field H_mesh] = ReadHDF5Dump([Sim_Path '/Hf.h5'], 'Range',{[],[],0});
% calc Bx,By, B1p, B1m normalize to the input-power
Bx = MUE0*H_field.FD.values{1}(:,:,1,1)/sqrt(P0_in);
By = MUE0*H_field.FD.values{1}(:,:,1,2)/sqrt(P0_in);
B1p = 0.5*(Bx+1j*By);
B1m = 0.5*(Bx-1j*By);
% create a 2D grid to plot on
[X Y] = ndgrid(H_mesh.lines{1},H_mesh.lines{2});

filter = sqrt(X.^2+Y.^2)<0.1;
Dump2VTK([Sim_Path '/B1p_xy.vtk'], abs(B1p).*filter, H_mesh, 'B-Field');
Dump2VTK([Sim_Path '/B1m_xy.vtk'], abs(B1m).*filter, H_mesh, 'B-Field');

% B1+ plot
figure()
h = pcolor(X,Y,log10(abs(B1p)));
set(h,'EdgeColor','none');
xlabel('x -->');
ylabel('y -->');
title('B_1^+ field (dB)');
axis equal tight

% B1- plot
figure()
h = pcolor(X,Y,log10(abs(B1m)));
set(h,'EdgeColor','none');
xlabel('x -->');
ylabel('y -->');
title('B_1^- field (dB)');
axis equal tight

%% dump to vtk to view in Paraview
Dump2VTK([Sim_Path '/SAR_xy.vtk'], SAR, SAR_mesh, 'SAR');

%%
ConvertHDF5_VTK([Sim_Path '/Hf.h5'],[Sim_Path '/B1_xy'], 'Range',{[],[],0}, 'weight', MUE0/sqrt(P0_in), 'FieldName', 'B1-field');

%% validation by calculating the power budget
% calculate the radiated power
nf2ff = CalcNF2FF(nf2ff, Sim_Path, f0, 0, 0);
p_rad  = nf2ff.Prad;

% read dissipated power in the phantom
p_loss = ReadHDF5Attribute([Sim_Path '/SAR.h5'],'/FieldData/FD/f0','power');

% display results, should add up to 100% (+/- 5% error margin)
disp([' power loss in phantom: ' num2str(p_loss) ' W (' num2str(p_loss/P0_in*100) '%)']);
disp([' power radiated: ' num2str(p_rad) ' W (' num2str(p_rad/P0_in*100) '%)']);
