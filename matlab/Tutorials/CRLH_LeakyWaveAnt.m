%
% Tutorials / CRLH Leaky Wave Antenna
%
% Tested with
%  - Octave 11.3
%  - openEMS v0.37
%
% (C) 2011-2026 Thorsten Liebig <thorsten.liebig@gmx.de>

close all
clear
clc

%% Simulation Parameters
%% ---------------------
%% Configure the CRLH leaky-wave antenna geometry and frequency sweep. The
%% CRLH struct fields (gap widths, stub dimensions, via radius) control the
%% left-handed capacitance and inductance that set the transition frequency;
%% N_Cells determines array length and therefore gain and beam-scanning range.
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

%% FDTD Parameters and Excitation
%% --------------------------------
%% A Gaussian pulse spanning f_start to f_stop excites all frequencies of
%% interest in a single run. PML absorbing boundaries on all six faces
%% model an open-space environment; EndCriteria stops the solver once the
%% stored energy has decayed to 0.1 % of its peak value.
FDTD = InitFDTD('EndCriteria', 1e-3);
FDTD = SetGaussExcite( FDTD, (f_start+f_stop)/2, (f_stop-f_start)/2 );
BC   = {'PML_8' 'PML_8' 'PML_8' 'PML_8' 'PML_8' 'PML_8'};
FDTD = SetBoundaryCond( FDTD, BC );

%% Mesh Setup and CRLH Unit Cell Instantiation
%% ---------------------------------------------
%% Seed the mesh with coarse lines at the domain boundaries and substrate
%% layer interfaces, then instantiate each of the N_Cells unit cells via
%% CreateCRLH, which appends fine lines around every metal edge. SmoothMesh
%% then interpolates a graded transition so the cell-size ratio stays below 1.5.
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

%% Substrate Layer Definition
%% ---------------------------
%% Three dielectric layers fill the region between the ground plane and the
%% top signal layer. Loss is modeled via an effective conductivity derived
%% from the loss tangent at 3 GHz; the layered stack realizes the
%% inter-layer coupling that produces the CRLH series and shunt reactances.
substratelines = [0 substratelines];
for n=1:numel(substrate_thickness)
    CSX = AddMaterial( CSX, ['substrate' int2str(n)] );
    CSX = SetMaterialProperty( CSX, ['substrate' int2str(n)], 'Epsilon', substrate_epsr(n), 'Kappa', substrate_tanD(n)*substrate_epsr(n)*EPS0*2*pi*3e9 );
    start = [-feed_length-(N_Cells*CRLH.LL)/2, -substrate_width/2, substratelines(n)];
    stop  = [+feed_length+(N_Cells*CRLH.LL)/2,  substrate_width/2, substratelines(n+1)];
    CSX = AddBox( CSX, ['substrate' int2str(n)], 0, start, stop );
end

%% Ground Plane and MSL Feed Ports
%% ---------------------------------
%% A full-width ground plane is placed at z = 0. Two 50-ohm microstrip
%% ports flank the CRLH array: port 1 is actively excited while port 2 is
%% matched-terminated, enabling both S11 (return loss) and S21 (insertion
%% loss) to be extracted from a single simulation run.
%ground plane
% priority 10, not 0: substrate1 spans z=0 as well, and on a priority tie the
% property added first wins -- at priority 0 the ground plane is shadowed by
% the substrate and never makes it into the operator ("Unused primitive").
CSX = AddMetal( CSX, 'ground' );
start = [-feed_length-(N_Cells*CRLH.LL)/2, -substrate_width/2, 0];
stop  = [+feed_length+(N_Cells*CRLH.LL)/2,  substrate_width/2, 0];
CSX = AddBox( CSX, 'ground', 10, start, stop );

CSX = AddMetal( CSX, 'PEC' );
portstart = [ -feed_length-(N_Cells*CRLH.LL)/2 , -CRLH.LW/2, substratelines(end)];
portstop  = [ -(N_Cells*CRLH.LL)/2,  CRLH.LW/2, 0];
[CSX,port{1}] = AddMSLPort( CSX, 999, 1, 'PEC', portstart, portstop, 0, [0 0 -1], 'ExcitePort', true, 'MeasPlaneShift',  feed_length/2, 'Feed_R', 50);

portstart = [ feed_length+(N_Cells*CRLH.LL)/2 , -CRLH.LW/2, substratelines(end)];
portstop  = [ +(N_Cells*CRLH.LL)/2,   CRLH.LW/2, 0];
[CSX,port{2}] = AddMSLPort( CSX, 999, 2, 'PEC', portstart, portstop, 0, [0 0 -1], 'MeasPlaneShift',  feed_length/2, 'Feed_R', 50 );

%% Near-Field to Far-Field Box
%% ----------------------------
%% The NF2FF box records tangential E and H field components on a closed
%% Huygens surface surrounding the antenna. It must sit at least 10 cells
%% inside the PML to avoid corrupted samples; OptResolution controls the
%% surface sampling density and therefore the far-field integration accuracy.
start = [mesh.x(1)   mesh.y(1)   mesh.z(1)  ] + 10*resolution;
stop  = [mesh.x(end) mesh.y(end) mesh.z(end)] - 10*resolution;
[CSX nf2ff] = CreateNF2FFBox(CSX, 'nf2ff', start, stop, 'OptResolution', nf2ff_resolution);

%% Write Geometry, Preview, and Run Simulation
%% ---------------------------------------------
%% Serialize the complete CSXCAD structure to an XML file that openEMS reads
%% at startup. CSXGeomPlot launches AppCSXCAD so the geometry can be
%% visually verified before committing to the approximately 30-minute run.
Sim_Path = 'tmp_CRLH_LeakyWave';
Sim_CSX = 'CRLH.xml';

CleanupSimPath(Sim_Path);

WriteOpenEMS( [Sim_Path '/' Sim_CSX], FDTD, CSX );
CSXGeomPlot( [Sim_Path '/' Sim_CSX] );
RunOpenEMS( Sim_Path, Sim_CSX );

%% S-Parameter Post-Processing
%% ----------------------------
%% calcPort reads time-domain probe data and applies a DFT over the
%% specified frequency vector; RefPlaneShift de-embeds the feed lines so
%% the reference plane sits at the CRLH array boundary. S11 well below
%% -10 dB across 1.9-4.2 GHz confirms good impedance matching of the LWA.
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

%% 3D Radiation Pattern Calculation
%% ----------------------------------
%% CalcNF2FF transforms the recorded Huygens-surface fields into far-field
%% E-theta and E-phi components at each frequency in f_rad. The sweep from
%% 1.9 to 4.2 GHz captures the frequency-dependent beam-scanning behavior
%% that is the defining characteristic of a CRLH leaky-wave antenna.
phi = 0:2:360;
theta = 0:2:180;

disp( 'calculating 3D far field pattern...' );
nf2ff = CalcNF2FF(nf2ff, Sim_Path, f_rad, theta*pi/180, phi*pi/180, 'Outfile','3D_Pattern.h5', 'Mode', 0,'Verbose',1);

%% Directivity and Radiation Efficiency
%% --------------------------------------
%% Plot maximum directivity (dBi) and radiation efficiency on dual y-axes
%% over the beam-scanning frequency range. Efficiency is the ratio of
%% radiated power Prad to accepted input power P_in, revealing how much
%% energy is dissipated in the substrate versus radiated into free space.
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

%% VTK Far-Field Pattern Export
%% ------------------------------
%% Write one VTK file per frequency for 3D visualization in ParaView.
%% Normalizing each pattern by its own peak and scaling by Dmax produces a
%% directivity-weighted surface, making beam-scanning and sidelobe evolution
%% across the frequency sweep immediately visible when loaded as an animation.
disp( 'dumping 3D far field pattern to vtk, use Paraview to visualize...' );
for n=1:numel(f_rad)
    E_far_normalized_3D = nf2ff.E_norm{n} / max(max(nf2ff.E_norm{n})) * nf2ff.Dmax(n);
    DumpFF2VTK( [Sim_Path '/FF_Pattern_' int2str(f_rad(n)/1e6) 'MHz.vtk'],E_far_normalized_3D,theta,phi,'scale',1e-3);
end

