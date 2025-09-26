close all
clear
clc

%confirm_recursive_rmdir(0);  % Uncomment to enable system removing directory without asking

% Frequency
f0 = 500.13e6;
fc = 200e6;   % Ignored when exciteMode is sinus

%% Setup the simulation
physical_constants; %get some physical constants like c0 and MUE0
unit = 1e-3; % all length in mm
max_res = c0 / (f0+fc) / unit / 20;
Airbox = c0 / (f0-fc) / unit / 25;

%% Saddle-Coil parameters
saddle.height = 20;
saddle.diameter = 4;
saddle.radius = saddle.diameter / 2;
saddle.angle = 90 * pi / 180;
saddle.wireRadius = 0.25;
saddle.wireSpace = 0.25;
saddle.wireSpaceFromCenter = 2*saddle.wireRadius + saddle.wireSpace;
saddle.rotation = saddle.angle / 2;   % Rotation along Z-axis
saddle.material = 1;  % 0: Perfect Electrical Conductor (PEC)
                      % 1: Copper

% Excitation port parameters
port_resist = 1000;

% Enable parameters
enableMatchingCapa = 0;
enableTuningCapa = 0;
enableBoreShield = 1;
enableGeometricPlot = 1;
enableStartSimulation = 0;
exciteMode = 1;     % 0: Sinus (f = f0)
                    % 1: Gaussian (f0-fc < f < f0+fc)

% Tuning Matching parameters
capaM_value = 16e-12;
capaT_value = 1.2e-12;

% Bore parameters
bore.radius = 10;
bore.shieldThickness = 0.5;

% Mesh resolution
mesh_res = [max_res max_res max_res];


%%% FDTD and CSX initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exciteMode
  FDTD = InitFDTD('NrTS',400e3,'EndCriteria',1e-5); %
  FDTD = SetGaussExcite(FDTD, f0, fc);
else
  FDTD = InitFDTD('NrTS',15e3,'EndCriteria',1e-4);
  FDTD = SetSinusExcite(FDTD, f0);
endif

% boundary conditions
BC = {'MUR' 'MUR' 'MUR' 'MUR' 'PML_8' 'PML_8'}; %pml in pos. and neg. z-direction
FDTD = SetBoundaryCond(FDTD,BC);
CSX = InitCSX();


%%% Building the system to simulate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
angularWireSpace = saddle.wireSpace / saddle.radius;
angularWireSpaceFromCenter = saddle.wireSpaceFromCenter / saddle.radius;


%% Building the saddle
a = (pi+saddle.angle-angularWireSpaceFromCenter)/2 - saddle.rotation;
points(:,1) = [cos(a)*(saddle.radius+saddle.wireSpaceFromCenter) sin(a)*(saddle.radius+saddle.wireSpaceFromCenter)+saddle.wireSpaceFromCenter 0];
count = 1;
angles = linspace((pi+saddle.angle-angularWireSpaceFromCenter)/2, saddle.angle-angularWireSpaceFromCenter, 9);
for a = angles
  count = count + 1;
  a = a - saddle.rotation;
  points(:,count) = [cos(a)*(saddle.radius+saddle.wireSpaceFromCenter) sin(a)*(saddle.radius+saddle.wireSpaceFromCenter) 0];
end
angles = linspace(saddle.angle-angularWireSpaceFromCenter, 0, 9);
for a = angles
  count = count + 1;
  a = a - saddle.rotation;
  points(:,count) = [cos(a)*saddle.radius sin(a)*saddle.radius 0];
end
angles = linspace(0, saddle.angle, 9);
for a = angles
  count = count + 1;
  a = a - saddle.rotation;
  points(:,count) = [cos(a)*saddle.radius sin(a)*saddle.radius saddle.height];
end
angles = linspace(saddle.angle, pi+saddle.angle, 13);
for a = angles
  count = count + 1;
  a = a - saddle.rotation;
  points(:,count) = [cos(a)*saddle.radius sin(a)*saddle.radius 0];
end
angles = linspace(pi+saddle.angle, pi, 9);
for a = angles
  count = count + 1;
  a = a - saddle.rotation;
  points(:,count) = [cos(a)*saddle.radius sin(a)*saddle.radius saddle.height];
end
count = count + 1;
points(:,count) = [cos(pi - saddle.rotation)*saddle.radius sin(pi - saddle.rotation)*saddle.radius saddle.wireSpaceFromCenter];
count = count + 1;
points(:,count) = [cos(pi - saddle.rotation)*(saddle.radius+saddle.wireSpaceFromCenter) sin(pi - saddle.rotation)*(saddle.radius+saddle.wireSpaceFromCenter) saddle.wireSpaceFromCenter];
angles = linspace(pi, (pi+saddle.angle+angularWireSpaceFromCenter)/2, 9);
for a = angles
  count = count + 1;
  a = a - saddle.rotation;
  points(:,count) = [cos(a)*(saddle.radius+saddle.wireSpaceFromCenter) sin(a)*(saddle.radius+saddle.wireSpaceFromCenter) 0];
end
points(:,count+1) = [points(:,end)(1) points(:,end)(2)+saddle.wireSpaceFromCenter points(:,end)(3)];

points = points.*1000;  % To correct imprecision of cos and sin functions
points = round(points);
points = points./1000;

if saddle.material
  CSX = AddMaterial(CSX,'saddle');
  CSX = SetMaterialProperty(CSX,'saddle','Kappa',56e6);
else
  CSX = AddMetal(CSX, 'saddle');
endif
CSX = AddWire(CSX, 'saddle', 10, points, saddle.wireRadius);


%% Adding the capacitors and excitation port
if enableTuningCapa
  % Add tunning capacitor
  start = points(:,2);
  stop = points(:,end-1);
  CSX = AddLumpedElement(CSX, 'CapaT', 0, 'Caps', 1, 'C', capaT_value);
  CSX = AddBox(CSX, 'CapaT', 0, [start(1) start(2)-0.25 start(3)-0.25], [stop(1) stop(2)+0.25 stop(3)+0.25]);
endif

if enableMatchingCapa
  % Add matching capacitor
  capaM_length = 0.5;
  capaM_start = [points(:,end)(1)-saddle.wireRadius points(:,end)(2) points(:,end)(3)-saddle.wireRadius];
  capaM_stop = [points(:,end)(1)+saddle.wireRadius points(:,end)(2)+saddle.wireSpaceFromCenter points(:,end)(3)+saddle.wireRadius];
  CSX = AddLumpedElement(CSX, 'CapaM', 1, 'Caps', 1, 'C', capaM_value);
  CSX = AddBox(CSX, 'CapaM', 0, capaM_start, capaM_stop);

  % Extending saddle start for capaM
  points_capaM_extend(:,1) = points(:,1);
  points_capaM_extend(:,2) = [points(:,1)(1) points(:,1)(2)+saddle.wireSpaceFromCenter points(:,1)(3)];
  CSX = AddWire(CSX, 'saddle', 10, points_capaM_extend, saddle.wireRadius);

  % Set excitation port start/stop
  start = points_capaM_extend(:,2);
  stop = [capaM_stop(1)-saddle.wireRadius capaM_stop(2) capaM_stop(3)-saddle.wireRadius];
else
  % Set excitation port start/stop
  start = points(:,1);
  stop = points(:,end);
endif

% Add excitation port
[CSX port] = AddLumpedPort(CSX, 100, 1, port_resist, start, stop, [1 0 0], true);


%% Add independent voltage and current probes
% Voltage probe
CSX = AddProbe(CSX, 'ut1', 0);
CSX = AddBox(CSX, 'ut1', 0, stop, start);

% Current probe
refPoint = start;
CSX = AddProbe(CSX, 'it1', 1);
start = [refPoint(1)-(saddle.wireRadius+saddle.wireSpace/2) refPoint(2)-saddle.wireSpace refPoint(3)-(saddle.wireRadius+saddle.wireSpace/2)];
stop = [refPoint(1)+(saddle.wireRadius+saddle.wireSpace/2) refPoint(2)-saddle.wireSpace refPoint(3)+(saddle.wireRadius+saddle.wireSpace/2)];
CSX = AddBox(CSX,'it1', 0 ,start,stop);


%% Add bore shield
if enableBoreShield
  start = [0 0 -saddle.wireRadius-Airbox];
  stop = [0 0 saddle.height+saddle.wireRadius+Airbox];
  CSX = AddMetal(CSX, 'BoreShield');
  CSX = AddCylindricalShell(CSX, 'BoreShield', 10, start, stop, bore.radius, bore.shieldThickness);
end


%%% Building the mesh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mesh = DetectEdges(CSX);

mesh.x = [mesh.x -saddle.radius-saddle.wireRadius saddle.radius+saddle.wireRadius];
mesh.y = [mesh.y -saddle.radius-saddle.wireRadius saddle.radius+saddle.wireRadius];
mesh.z = [mesh.z -saddle.wireRadius saddle.height+saddle.wireRadius];

mesh.x = sort(mesh.x);
mesh.y = sort(mesh.y);
mesh.z = sort(mesh.z);

tempMesh = mesh;
mesh.x = [];
mesh.y = [];
mesh.z = [];

for i = 1:length(tempMesh.x)-1
  N = round((tempMesh.x(i+1) - tempMesh.x(i)) / saddle.wireRadius);
  mesh.x = [mesh.x linspace(tempMesh.x(i), tempMesh.x(i+1), N)];
end

for i = 1:length(tempMesh.y)-1
  N = round((tempMesh.y(i+1) - tempMesh.y(i)) / saddle.wireRadius);
  mesh.y = [mesh.y linspace(tempMesh.y(i), tempMesh.y(i+1), N)];
end

for i = 1:length(tempMesh.z)-1
  N = round((tempMesh.z(i+1) - tempMesh.z(i)) / saddle.wireRadius);
  mesh.z = [mesh.z linspace(tempMesh.z(i), tempMesh.z(i+1), N)];
end

mesh.x = [min(mesh.x)-Airbox mesh.x max(mesh.x)+Airbox];
mesh.y = [min(mesh.y)-Airbox mesh.y max(mesh.y)+Airbox];
mesh.z = [min(mesh.z)-Airbox mesh.z max(mesh.z)+Airbox];

mesh = SmoothMesh(mesh, mesh_res, 1.4);
mesh = AddPML(mesh, [0 0 0 0 8 8]);

CSX = DefineRectGrid(CSX, unit, mesh);



%%% Add field analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
start = [-2*saddle.radius -2*saddle.radius saddle.height/2];
stop  = [2*saddle.radius 2*saddle.radius saddle.height/2];
CSX = AddDump(CSX,'Ht_xy','DumpType',1,'FileType',1);
CSX = AddBox(CSX,'Ht_xy',0, start, stop);

start = [0 -2*saddle.radius 0];
stop  = [0 2*saddle.radius saddle.height];
CSX = AddDump(CSX,'Ht_zy','DumpType',1,'FileType',1);
CSX = AddBox(CSX,'Ht_zy',0, start, stop);


%%% Setup files for recording simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sim_Path = 'tmp';
Sim_CSX = 'saddle.xml';

[status, message, messageid] = rmdir(Sim_Path, 's'); % clear previous directory
[status, message, messageid] = mkdir(Sim_Path ); % create empty simulation folder

WriteOpenEMS([Sim_Path '/' Sim_CSX], FDTD, CSX);
if enableGeometricPlot
  CSXGeomPlot([Sim_Path '/' Sim_CSX]);  % Open 3D viewer to show geometry
endif
if enableStartSimulation
  RunOpenEMS(Sim_Path, Sim_CSX);        % Start simulation
endif


%%% POST-PROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exciteMode
  f = linspace(f0-fc, f0+fc, 501);
  port = calcPort(port, Sim_Path, f);
else
  f = f0;
  port = calcPort(port, Sim_Path, f, 'SignalType', 'periodic');
endif

U = ReadUI('ut1','tmp/');
I = ReadUI('it1','tmp/');

if exciteMode
  Zin = port.uf.tot ./ port.if.tot;
  L = imag(Zin)./(f*2*pi);
  R = real(Zin);
  s11 = port.uf.ref ./ port.uf.inc;

  subplot(2,1,1);
  plot(f*1e-6,L*1e9,'Linewidth',2);
  xlabel('Frequency (MHz)');
  ylabel('Coil inductance (nH)');
  grid on;
  subplot(2,1,2);
  plot(f*1e-6,R,'Linewidth',2);
  hold on
  plot(f*1e-6,imag(Zin),'r','Linewidth',2);
  xlabel('Frequency (MHz)');
  ylabel('Resistance (Ohm)');
  grid on;
  legend( {'real','imaginary'}, 'location', 'northwest' )

  figure
  plot( f/1e6, 20*log10(abs(s11)), 'k-', 'Linewidth', 2 );
  ylim([-40 10]);
  grid on
  title( 'Reflection coefficient S_{11}' );
  ylabel( 'Reflection coefficient |S_{11}|' );
  xlabel( 'Frequency (MHz)' );
endif

figure
subplot(2,1,1);
plot(port.ut.time/1e-6,port.ut.tot,'Linewidth',2);
hold on
plot(U.TD{1}.t/1e-6, U.TD{1}.val, 'Linewidth',2);
xlabel('Time (us)');
ylabel('Amplitude (V)');
grid on;
legend('Port', 'Independent probe')
subplot(2,1,2);
plot(port.it.time/1e-6,port.it.tot,'Linewidth',2);
hold on
plot(I.TD{1}.t/1e-6, I.TD{1}.val, 'Linewidth',2);
xlabel('Time (us)');
ylabel('Amplitude (A)');
grid on;
legend('Port', 'Independent probe')

if exciteMode
  Freq = input('Enter resonant frequency : ');
  if ~ismember(Freq, f)
    delta = abs(Freq - f);
    closestFreq = f(find(delta == min(delta)));
    disp([int2str(Freq) ' is not in the frequency list. Taking the nearest value : ' int2str(closestFreq(1))]);
    Freq = closestFreq(1);
  endif
else
  Freq = f;
endif

disp(['Dumping resonant H-field XY at f = ' int2str(Freq) ' to vtk file, use Paraview to visualize']);
ConvertHDF5_VTK([Sim_Path '/Ht_xy.h5'],[Sim_Path '/Hf_xy_'],'Frequency',Freq,'FieldName','H-Field');
disp(['Dumping resonant H-field ZY at f = ' int2str(Freq) ' to vtk file, use Paraview to visualize']);
ConvertHDF5_VTK([Sim_Path '/Ht_zy.h5'],[Sim_Path '/Hf_zy_'],'Frequency',Freq,'FieldName','H-Field');
