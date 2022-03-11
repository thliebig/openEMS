%
% Tutorials / CPW_Line
%
% Description at:
%
% Tested with
%  - Octave 3.8.1
%  - openEMS v0.0.32
%
% (C) 2014 Thorsten Liebig <thorsten.liebig@gmx.de>

close all
clear
clc

%% setup the simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
physical_constants;
unit = 1e-6; % specify everything in um
CPW_length = 40000;
CPW_port_length = 10000;
CPW_width  = 1000;
CPW_gap    = 140;
substrate_thickness = 512;
substrate_width = 5000
substrate_epr = 3.66;
f_max = 10e9;
air_spacing = 7000

% use a finite line CPW waveguide
if 1
  feed_R = 50;
  pml_add_cells = [8 8 8 8 8 8];
  feed_shift_cells = 0;
  x_spacing = air_spacing;
else % or use a waveguide with start/end in a pml
  feed_R = inf; % CPW ends in a pml --> disable termination resistance 
  feed_shift_cells = 10; % CPW ends in an 8 cells thick pml --> shift feed 10 cells
  pml_add_cells = [0 0 8 8 8 8]; % do not add air-space in x-direction
  x_spacing = 0; % do not add air-space in x-direction
end

%% setup FDTD parameters & excitation function %%%%%%%%%%%%%%%%%%%%%%%%%%%%
FDTD = InitFDTD('EndCriteria', 1e-4);
FDTD = SetGaussExcite( FDTD, f_max/2, f_max/2 );
BC   = [2 2 2 2 2 2];
FDTD = SetBoundaryCond( FDTD, BC );

%% setup CSXCAD geometry & mesh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CSX = InitCSX();
resolution = c0/(f_max*sqrt(substrate_epr))/unit /30; % resolution of lambda/50
edge_res = 40;
mesh.x = SmoothMeshLines( [0 CPW_length/2 CPW_length/2+x_spacing], resolution, 1.5 ,0 );
mesh.x = unique(sort([-mesh.x mesh.x]));
mesh.y = SmoothMeshLines( [CPW_width/2+[-edge_res/3 +edge_res/3*2] CPW_gap+CPW_width/2+[-edge_res/3*2 +edge_res/3]], edge_res , 1.5 ,0);
mesh.y = SmoothMeshLines( [0 mesh.y], edge_res*2, 1.3 ,0);
mesh.y = SmoothMeshLines( [0 mesh.y substrate_width/2 substrate_width/2+air_spacing], resolution, 1.3 ,0);
mesh.y = unique(sort([-mesh.y mesh.y]));
mesh.z = SmoothMeshLines( [-air_spacing linspace(0,substrate_thickness,5) substrate_thickness+air_spacing], resolution );

mesh = AddPML(mesh, pml_add_cells);
CSX = DefineRectGrid( CSX, unit, mesh );

%% substrate
CSX = AddMaterial( CSX, 'RO4350B' );
CSX = SetMaterialProperty( CSX, 'RO4350B', 'Epsilon', substrate_epr );
start = [-CPW_length/2, -substrate_width/2, 0];
stop  = [+CPW_length/2, +substrate_width/2, substrate_thickness];
CSX = AddBox( CSX, 'RO4350B', 0, start, stop );

%%
CSX = AddMetal( CSX, 'CPW_PORT' );

%% CPW port, with the measurement plane at the end of each port
portstart = [ -CPW_length/2                , -CPW_width/2, substrate_thickness];
portstop  = [ -CPW_length/2+CPW_port_length,  CPW_width/2, substrate_thickness];
[CSX,port{1}] = AddCPWPort( CSX, 999, 1, 'CPW_PORT', portstart, portstop, CPW_gap, 'x', [0 1 0], 'ExcitePort', true, 'FeedShift', feed_shift_cells*resolution, 'MeasPlaneShift',  CPW_port_length, 'Feed_R', feed_R);

portstart = [ CPW_length/2                , -CPW_width/2, substrate_thickness];
portstop  = [ CPW_length/2-CPW_port_length,  CPW_width/2, substrate_thickness];
[CSX,port{2}] = AddCPWPort( CSX, 999, 2, 'CPW_PORT', portstart, portstop, CPW_gap, 'x', [0 1 0], 'MeasPlaneShift',  CPW_port_length, 'Feed_R', feed_R);

%% CPW
CSX = AddMetal( CSX, 'CPW');
start = [ -CPW_length/2+CPW_port_length, -CPW_width/2, substrate_thickness];
stop  = [ +CPW_length/2-CPW_port_length,  CPW_width/2, substrate_thickness];
CSX = AddBox(CSX, 'CPW', 999, start, stop);

%% CPW grounds
CSX = AddMetal( CSX, 'GND' );
start = [-CPW_length/2, -CPW_width/2-CPW_gap, substrate_thickness];
stop  = [+CPW_length/2, -substrate_width/2  , substrate_thickness];
CSX = AddBox(CSX, 'GND', 999, start, stop);

start = [-CPW_length/2, +CPW_width/2+CPW_gap, substrate_thickness];
stop  = [+CPW_length/2, +substrate_width/2  , substrate_thickness];
CSX = AddBox(CSX, 'GND', 999, start, stop);

%% write/show/run the openEMS compatible xml-file
Sim_Path = 'tmp';
Sim_CSX = 'CPW.xml';

[status, message, messageid] = rmdir( Sim_Path, 's' ); % clear previous directory
[status, message, messageid] = mkdir( Sim_Path ); % create empty simulation folder

WriteOpenEMS( [Sim_Path '/' Sim_CSX], FDTD, CSX );
CSXGeomPlot( [Sim_Path '/' Sim_CSX] );
RunOpenEMS( Sim_Path, Sim_CSX );

%% post-processing
close all
f = linspace( 1e6, f_max, 1601 );
port = calcPort( port, Sim_Path, f, 'RefImpedance', 50);

s11 = port{1}.uf.ref./ port{1}.uf.inc;
s21 = port{2}.uf.ref./ port{1}.uf.inc;

plot(f/1e9,20*log10(abs(s11)),'k-','LineWidth',2);
hold on;
grid on;
plot(f/1e9,20*log10(abs(s21)),'r--','LineWidth',2);
legend('S_{11}','S_{21}');
ylabel('S-Parameter (dB)','FontSize',12);
xlabel('frequency (GHz) \rightarrow','FontSize',12);


