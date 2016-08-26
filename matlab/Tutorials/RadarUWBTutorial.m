% Tutorial on time delay and signal integrity for radar 
% and UWB applications
% 
% Author: Georg Michel, 2016

 clear;
 close all;

physical_constants;

% --- start of configuration section ---

% In radar and ultrawideband applications it is important to know the
% delay and fidelity of RF pulses. The delay is the retardation of the
% signal from the source to the phase center of the antenna. It is
% composed out of linear delay, dispersion and minimum-phase
% delay. Dispersion due to waveguides or frequency-dependent
% permittivity and minimum-phase delay due to resonances will degrade
% the fidelity which is the normalized similarity between excitation and
% radiated signal. In this tutorial you can examine the performance of a
% simple ultrawideband (UWB) monopole. The delay and fidelity of this
% antenna are calculated and plotted. You can compare these properties
% in different channels.
% 
% The Gaussian excitation is set to the same 3dB bandwidth as the
% channels of the IEEE 802.15.4 UWB PHY. One exeption is channel4twice
% which has the double bandwidth of channel 4. It can be seen that the
% delay is larger and the fidelity is smaller in the vicinity of the
% (undesired) resonances of the antenna. Note that for a real UWB system
% the total delay and fidelity result from both the transmitting and
% receiving antenna or twice the delay and the square of the fidelity
% for monostatic radar. 
% 
% The resolution of the delay will depend on the 'Oversampling'
% parameter to InitFDTD.  See the description of DelayFidelity
% 
% In the configuration section below you can uncomment the respective
% parameter settings. 


% 3dB bandwidth is 0.3925 times 20dB bandwidth for Gaussian excitation

%suffix = "channel1";
%f_0 = 3.5e9; % center frequency of the channel
%f_c = 0.25e9 / 0.3925; 

%suffix = "channel2";
%f_0 = 4.0e9; % center frequency of the channel
%f_c = 0.25e9 / 0.3925; 

%suffix = "channel4";
%f_0 = 4.0e9; % center frequency of the channel
%f_c = 0.5e9 / 0.3925; 

%suffix = "channel5";
%f_0 = 6.5e9; % center frequency of the channel
%f_c = 0.25e9 / 0.3925; 

%suffix = "channel7";
%f_0 = 6.5e9; % center frequency of the channel
%f_c = 0.5e9 / 0.3925; 

suffix = "channel4twice";
f_0 = 4.0e9; % center frequency of the channel
f_c = 1e9 / 0.3925; 

tilt = 45 * pi / 180; % polarization tilt angle against co-polarization (90DEG is cross polarized) 

% --- end of configuration section ---


Sim_Path = 'tmp';
Sim_CSX = 'uwb.xml';


substrate.epsR = 4.2; % FR4
substrate.height = 0.707;



unit = 1e-3;


fineResolution = C0 / (f_0 + f_c) / sqrt(substrate.epsR) / unit / 120;
coarseResolution = C0/(f_0 + f_c) / unit / 60;
CSX = InitCSX();
CSX = AddMetal(CSX, 'Ground' );
CSX = AddMetal(CSX, 'Patch');
CSX = AddMaterial(CSX, 'Substrate' );
CSX = SetMaterialProperty(CSX, 'Substrate', 'Epsilon', substrate.epsR);



CSX = AddBox(CSX, 'Ground', 2, [-16, -16, -substrate.height], [16, 0, -substrate.height]);
CSX = AddBox(CSX, 'Patch', 2, [-1.15, -16, 0], [1.15, 8.5, 0]);
CSX = AddBox(CSX, 'Patch', 2, [-7, 0.5, 0], [7, 14.5, 0]);


mesh.x = [];
mesh.y = [];
mesh.z = [-substrate.height, 0];
mesh = DetectEdges(CSX, mesh, '2D_Metal_Edge_Res', fineResolution);
mesh.y(find(diff(mesh.y) < 0.6 * fineResolution) + 1) = []; % remove mesh lines which are too close
mesh.y = [mesh.y, 18]; % top of substrate
mesh = SmoothMesh(mesh, fineResolution);
mesh.x = [mesh.x -60, 60];
mesh.y = [mesh.y, -60, 65];
mesh.z = [mesh.z, -46, 45];
mesh = SmoothMesh(mesh, coarseResolution);
CSX = DefineRectGrid( CSX, unit, mesh);
[CSX, port] = AddLumpedPort( CSX, 999, 1, 50, [-1.15, -16, -substrate.height], [1.15, -15, 0], [0 0 1], true);
CSX = AddBox(CSX, 'Substrate', 1, [-16, -16, -substrate.height], [16, 18, 0]);
start = [mesh.x(10)     mesh.y(10)     mesh.z(10)];
stop  = [mesh.x(end-9) mesh.y(end-9) mesh.z(end-9)];
[CSX nf2ff] = CreateNF2FFBox(CSX, 'nf2ff', start, stop);

FDTD = InitFDTD( 'NrTs', 120000, 'EndCriteria', 1e-5, 'OverSampling', 20);
FDTD = SetGaussExcite(FDTD, f_0, f_c );
BC   = {'PML_8' 'PML_8' 'PML_8' 'PML_8' 'PML_8' 'PML_8'};
FDTD = SetBoundaryCond(FDTD, BC );


% remove old data, show structure, calculate new data
confirm_recursive_rmdir (0);
[status, message, messageid] = rmdir( Sim_Path, 's' ); % clear previous directory
[status, message, messageid] = mkdir( Sim_Path ); % create empty simulation folder
WriteOpenEMS([Sim_Path '/' Sim_CSX], FDTD, CSX);
CSXGeomPlot([Sim_Path '/' Sim_CSX]);
%return;
RunOpenEMS( Sim_Path, Sim_CSX);
%return;

freq  = linspace(f_0-f_c, f_0+f_c, 100);
port = calcPort(port, Sim_Path, freq);
s11 = port.uf.ref ./ port.uf.inc;

s11phase = unwrap(arg(s11));

%
% plot reflection coefficient S11
figure%("visible", "off");
ax = plotyy( freq/1e6, 20*log10(abs(s11)), freq/1e6, s11phase);
grid on
title( 'reflection coefficient S_{11}' );
xlabel( 'frequency f / MHz' );
ylabel( ax(1), 'reflection coefficient |S_{11}|' );
ylabel(ax(2), 'S_{11} phase (rad)');


phi = [0] * pi / 180;
theta = [-180:10:180] * pi / 180;
nf2ff = CalcNF2FF(nf2ff, Sim_Path, [f_0], theta, phi, 'Mode', 1);
figure%("visible", "off");
polarFF(nf2ff, 'xaxis', 'theta', 'freq_index', 1, 'logscale', -10);

[delay, fidelity] = DelayFidelity(nf2ff, port, Sim_Path, sin(tilt), cos(tilt), theta, phi, f_0, f_c, 'Mode', 1);

figure%("visible", "off");
plot(theta.'/pi*180, delay(:,1) * 3e11);
ylim([0, 80]);
xlabel('theta/DEG');
ylabel('delay/mm');
title("Delay Length");

figure%("visible", "off");
plot(theta.'/pi*180, fidelity(:,1));
ylim([0.95, 1]);
xlabel('theta/DEG');
ylabel('fidelity/a.u.');
title('Fidelity Factor');


print(1, ["s11_", suffix, ".png"]);
print(2, ["farfield_", suffix, ".png"]);
print(3, ["delay_mm_", suffix, ".png"]);
print(4, ["fidelity_", suffix, ".png"]);
return;

