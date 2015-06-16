%
% Tutorials / CRLH_LeakyWaveAnt
%
% Describtion at:
% http://openems.de/index.php/Tutorial:_CRLH_Leaky_Wave_Antenna
%
% Tested with
%  - Matlab 2011a / Octave 4.0
%  - openEMS v0.0.33
%
% (C) 2011-2015 Thorsten Liebig <thorsten.liebig@gmx.de>

close all
clear
clc

%% setup the simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
physical_constants;
unit = 1e-6; % specify everything in um

feed_length = 20000;

substrate_thickness = [1524 101 254];
substrate_epsr = [3.48 3.48 3.48];
substrate_tanD = [1 1 1]*1e-3;

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
Air_Spacer = 30000;

% frequency range of interest
f_start = 1e9;
f_stop  = 6e9;

% frequencies to calculate the 3D radiation pattern
f_rad = (1.9:0.05:4.2)*1e9;
nf2ff_resolution = c0/max(f_rad)/unit/15;

%% setup FDTD parameters & excitation function %%%%%%%%%%%%%%%%%%%%%%%%%%%%
FDTD = InitFDTD('EndCriteria', 1e-3);
FDTD = SetGaussExcite( FDTD, (f_start+f_stop)/2, (f_stop-f_start)/2 );
BC   = {'PML_8' 'PML_8' 'PML_8' 'PML_8' 'PML_8' 'PML_8'};
FDTD = SetBoundaryCond( FDTD, BC );

%% Setup a basic mesh and create the CRLH unit cell
CSX = InitCSX();
resolution = c0/(f_stop*sqrt(max(substrate_epsr)))/unit /30; % resolution of lambda/30

mesh.x = [-feed_length-(N_Cells*CRLH.LL)/2-Air_Spacer -feed_length-(N_Cells*CRLH.LL)/2 0 feed_length+(N_Cells*CRLH.LL)/2 feed_length+(N_Cells*CRLH.LL)/2+Air_Spacer];
mesh.y = [-Air_Spacer-substrate_width/2 0 Air_Spacer+substrate_width/2];
substratelines = cumsum(substrate_thickness);
mesh.z = [-0.5*Air_Spacer 0 cumsum(substrate_thickness) linspace(substratelines(end-1),substratelines(end),4) Air_Spacer];

% create the CRLH unit cells (will define additional fixed mesh lines)
pos_x = -(N_Cells*CRLH.LL)/2 + CRLH.LL/2;
for n=1:N_Cells
    [CSX mesh] = CreateCRLH(CSX, mesh, CRLH, resolution/4, [pos_x 0 0]);
    pos_x = pos_x + CRLH.LL;
end

% Smooth the given mesh
mesh = SmoothMesh(mesh, resolution, 1.5, 'algorithm',[1 3]);
CSX = DefineRectGrid( CSX, unit, mesh );

%% Setup the substrate layer
substratelines = [0 substratelines];
for n=1:numel(substrate_thickness)
    CSX = AddMaterial( CSX, ['substrate' int2str(n)] );
    CSX = SetMaterialProperty( CSX, ['substrate' int2str(n)], 'Epsilon', substrate_epsr(n), 'Kappa', substrate_tanD(n)*substrate_epsr(n)*EPS0*2*pi*3e9 );
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
[CSX,port{1}] = AddMSLPort( CSX, 999, 1, 'PEC', portstart, portstop, 0, [0 0 -1], 'ExcitePort', true, 'MeasPlaneShift',  feed_length/2, 'Feed_R', 50);

portstart = [ feed_length+(N_Cells*CRLH.LL)/2 , -CRLH.LW/2, substratelines(end)];
portstop  = [ +(N_Cells*CRLH.LL)/2,   CRLH.LW/2, 0];
[CSX,port{2}] = AddMSLPort( CSX, 999, 2, 'PEC', portstart, portstop, 0, [0 0 -1], 'MeasPlaneShift',  feed_length/2, 'Feed_R', 50 );

%% nf2ff calc
start = [mesh.x(1)   mesh.y(1)   mesh.z(1)  ] + 10*resolution;
stop  = [mesh.x(end) mesh.y(end) mesh.z(end)] - 10*resolution;
[CSX nf2ff] = CreateNF2FFBox(CSX, 'nf2ff', start, stop, 'OptResolution', nf2ff_resolution);

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
port = calcPort( port, Sim_Path, f, 'RefPlaneShift', feed_length*unit);

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

%% calculate 3D pattern
phi = 0:2:360;
theta = 0:2:180;

disp( 'calculating 3D far field pattern...' );
nf2ff = CalcNF2FF(nf2ff, Sim_Path, f_rad, theta*pi/180, phi*pi/180, 'Outfile','3D_Pattern.h5', 'Mode', 0,'Verbose',1);

%%
P_in = interp1(f, port{1}.P_acc, f_rad);

figure()

[AX,H1,H2] = plotyy(f_rad/1e9,nf2ff.Dmax',f_rad/1e9,100*nf2ff.Prad'./P_in,'plot');
grid on
xlabel( 'frequency (GHz)' );
set(get(AX(1),'Ylabel'),'String','directivity (dBi)')
set(get(AX(2),'Ylabel'),'String','radiation efficiency (%)')
set(H1,'Linewidth',2)
set(H2,'Linewidth',2)
set(H1,'Marker','*')
set(H2,'Marker','s')

drawnow

%%
disp( 'dumping 3D far field pattern to vtk, use Paraview to visualize...' );
for n=1:numel(f_rad)
    E_far_normalized_3D = nf2ff.E_norm{n} / max(max(nf2ff.E_norm{n})) * nf2ff.Dmax(n);
    DumpFF2VTK( [Sim_Path '/FF_Pattern_' int2str(f_rad(n)/1e6) 'MHz.vtk'],E_far_normalized_3D,theta,phi,'scale',1e-3);
end

