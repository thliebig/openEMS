%
% Tutorials / CRLH_LeakyWaveAnt
%
% Describtion at:
% http://openems.de/index.php/Tutorial:_CRLH_Leaky_Wave_Antenna
%
% Tested with
%  - Matlab 2011a / Octave 3.4.3
%  - openEMS v0.0.27
%
% (C) 2011,2012 Thorsten Liebig <thorsten.liebig@gmx.de>

close all
clear
clc

%% setup the simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
physical_constants;
unit = 1e-6; % specify everything in um

feed_length = 20000;

substrate_thickness = [1524 101 254];
substrate_epsr = [3.48 3.48 3.48];

N_Cells = 8;        %number of CRLH unit cells

CRLH.LL = 14e3;     %CRLH totel (line) length
CRLH.LW = 4e3;      %CRLH unit cell width (without the stubs)
CRLH.GLB = 1950;    %CRLH gap width bottom layer
CRLH.GLT = 4700;    %CRLH gap width top layer
CRLH.SL = 7800;     %CRLH stub length (bottom layer, both sides)
CRLH.SW = 1000;     %CRLH stub width  (bottom layer, both sides)
CRLH.VR = 250;      %CRLH via hole radius (stub -> ground)
CRLH.TopSig = sum(substrate_thickness);  %top layer height
CRLH.BottomSig = CRLH.TopSig - substrate_thickness(end);  %bottom layer height

substrate_width = CRLH.LW + 2*CRLH.SL;
Air_Spacer = 25000;

% frequency range of interest
f_start = 0.8e9;
f_stop  = 6e9;

f_rad = (1.9:0.1:4.2)*1e9;

Plot_3D_Rad_Pattern = 1; %this may take a long time! > 30min

%% setup FDTD parameters & excitation function %%%%%%%%%%%%%%%%%%%%%%%%%%%%
FDTD = InitFDTD( 20000 );
FDTD = SetGaussExcite( FDTD, (f_start+f_stop)/2, (f_stop-f_start)/2 );
BC   = {'PML_8' 'PML_8' 'PML_8' 'PML_8' 'PML_8' 'PML_8'};
FDTD = SetBoundaryCond( FDTD, BC );

%% Setup a basic mesh and create the CRLH unit cell
CSX = InitCSX();
resolution = c0/(f_stop*sqrt(max(substrate_epsr)))/unit /30; % resolution of lambda/30

mesh.x = [-feed_length-(N_Cells*CRLH.LL)/2-Air_Spacer -feed_length-(N_Cells*CRLH.LL)/2 0 feed_length+(N_Cells*CRLH.LL)/2 feed_length+(N_Cells*CRLH.LL)/2+Air_Spacer];
mesh.y = [-Air_Spacer-substrate_width/2 0 Air_Spacer+substrate_width/2];
substratelines = cumsum(substrate_thickness);
mesh.z = [-0.7*Air_Spacer 0 cumsum(substrate_thickness) linspace(substratelines(end-1),substratelines(end),4) Air_Spacer];

% create the CRLH unit cells (will define additional fixed mesh lines)
pos_x = -(N_Cells*CRLH.LL)/2 + CRLH.LL/2;
for n=1:N_Cells
    [CSX mesh] = CreateCRLH(CSX, mesh, CRLH, resolution/4, [pos_x 0 0]);
    pos_x = pos_x + CRLH.LL;
end

% Smooth the given mesh
mesh.x = SmoothMeshLines(mesh.x, resolution, 1.5, 0);
mesh.y = SmoothMeshLines(mesh.y, resolution, 1.5, 0);
mesh.z = SmoothMeshLines(mesh.z, resolution, 1.5, 0);
CSX = DefineRectGrid( CSX, unit, mesh );

%% Setup the substrate layer
substratelines = [0 substratelines];
for n=1:numel(substrate_thickness)
    CSX = AddMaterial( CSX, ['substrate' int2str(n)] );
    CSX = SetMaterialProperty( CSX, ['substrate' int2str(n)], 'Epsilon', substrate_epsr(n) );
    start = [-feed_length-(N_Cells*CRLH.LL)/2, -substrate_width/2, substratelines(n)];
    stop  = [+feed_length+(N_Cells*CRLH.LL)/2,  substrate_width/2, substratelines(n+1)];
    CSX = AddBox( CSX, ['substrate' int2str(n)], 0, start, stop );
end

%% add the feeding MSL ports
%ground plane
CSX = AddMetal( CSX, 'ground' );
start = [-feed_length-(N_Cells*CRLH.LL)/2, -substrate_width/2, 0];
stop  = [+feed_length+(N_Cells*CRLH.LL)/2,  substrate_width/2, 0];
CSX = AddBox( CSX, 'ground', 0, start, stop );

CSX = AddMetal( CSX, 'PEC' );
portstart = [ -feed_length-(N_Cells*CRLH.LL)/2 , -CRLH.LW/2, substratelines(end)];
portstop  = [ -(N_Cells*CRLH.LL)/2,  CRLH.LW/2, 0];
[CSX,portstruct{1}] = AddMSLPort( CSX, 999, 1, 'PEC', portstart, portstop, 0, [0 0 -1], 'ExcitePort', 'excite', 'MeasPlaneShift',  feed_length/2, 'Feed_R', 50);

portstart = [ feed_length+(N_Cells*CRLH.LL)/2 , -CRLH.LW/2, substratelines(end)];
portstop  = [ +(N_Cells*CRLH.LL)/2,   CRLH.LW/2, 0];
[CSX,portstruct{2}] = AddMSLPort( CSX, 999, 2, 'PEC', portstart, portstop, 0, [0 0 -1], 'MeasPlaneShift',  feed_length/2, 'Feed_R', 50 );


%% nf2ff calc
start = [mesh.x(1)   mesh.y(1)   mesh.z(1)  ] + 10*resolution;
stop  = [mesh.x(end) mesh.y(end) mesh.z(end)] - 10*resolution;
[CSX nf2ff] = CreateNF2FFBox(CSX, 'nf2ff', start, stop);

%% write/show/run the openEMS compatible xml-file
Sim_Path = 'tmp_CRLH_LeakyWave';
Sim_CSX = 'CRLH.xml';

[status, message, messageid] = rmdir( Sim_Path, 's' ); % clear previous directory
[status, message, messageid] = mkdir( Sim_Path ); % create empty simulation folder

WriteOpenEMS( [Sim_Path '/' Sim_CSX], FDTD, CSX );
CSXGeomPlot( [Sim_Path '/' Sim_CSX] );
RunOpenEMS( Sim_Path, Sim_CSX );

%% post-processing
close all
f = linspace( f_start, f_stop, 1601 );
port{1} = calcPort( portstruct{1}, Sim_Path, f, 'RefPlaneShift', feed_length*unit);
port{2} = calcPort( portstruct{2}, Sim_Path, f, 'RefPlaneShift', feed_length*unit);

s11 = port{1}.uf.ref./ port{1}.uf.inc;
s21 = port{2}.uf.ref./ port{1}.uf.inc;

plot(f/1e9,20*log10(abs(s11)),'k-','LineWidth',2);
hold on;
grid on;
plot(f/1e9,20*log10(abs(s21)),'r--','LineWidth',2);
l = legend('S_{11}','S_{21}','Location','Best');
set(l,'FontSize',12);
ylabel('S-Parameter (dB)','FontSize',12);
xlabel('frequency (GHz) \rightarrow','FontSize',12);
ylim([-40 2]);

drawnow

%% NFFF contour plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta = (0:3:359) - 180;
phi = [0 90];

disp( 'calculating far field at phi=[0 90] deg...' );

nf2ff = CalcNF2FF(nf2ff, Sim_Path, f_rad, theta*pi/180, phi*pi/180, 1, 'Verbose',1);
nf2ff = ReadNF2FF(nf2ff);

%%
% prepare figures
figure(10)
hold on;
grid on;
xlabel( 'theta (deg)' );
ylabel( 'directivity (dBi)');
title('phi = 0°');
ylim([-20 10]);
figure(11)
hold on;
grid on;
xlabel( 'theta (deg)' );
ylabel( 'directivity (dBi)');
title('phi = 90°');
ylim([-20 10]);
line_styles = {'b-','g:','r-.','c--','m-','y:','k-.'};

for n=1:numel(f_rad)
    f_res = f_rad(n)

    % display power and directivity
    disp( ['frequency: f = ' num2str(f_res/1e9) ' GHz']);
    disp( ['radiated power: Prad = ' num2str(nf2ff.Prad(n)) ' Watt']);
    disp( ['directivity: Dmax = ' num2str(nf2ff.Dmax(n)) ' (' num2str(10*log10(nf2ff.Dmax(n))) ' dBi)'] );

    % normalized directivity
    D_log = 20*log10(nf2ff.E_norm{n}/max(max(nf2ff.E_norm{n})));
    % directivity
    D_log = D_log + 10*log10(nf2ff.Dmax(n));

    figure(10)
    plot( nf2ff.theta, D_log(:,1) ,line_styles{1+mod(n-1,numel(line_styles))});
    hold on;

    figure(11)
    plot( nf2ff.theta, D_log(:,2) ,line_styles{1+mod(n-1,numel(line_styles))} );
    hold on;
end

if (Plot_3D_Rad_Pattern==0)
    return
end

%% calculate 3D pattern
phi = 0:3:360;
theta = 0:3:180;

disp( 'calculating 3D far field pattern...' );
nf2ff = CalcNF2FF(nf2ff, Sim_Path, f_rad, theta*pi/180, phi*pi/180, 1, 'Verbose',2);
nf2ff = ReadNF2FF(nf2ff);

%%
disp( 'dumping 3D far field pattern to vtk, use Paraview to visualize...' );
for n=1:numel(f_rad)
    E_far_normalized_3D = nf2ff.E_norm{n} / max(max(nf2ff.E_norm{n})) * nf2ff.Dmax(n);
    DumpFF2VTK( [Sim_Path '/FF_Pattern_' int2str(f_rad(n)/1e6) 'MHz.vtk'],E_far_normalized_3D,theta,phi,1e-3);
end

