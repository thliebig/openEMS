%
% Tutorials / Rect_Waveguide
%
% Describtion at:
% http://openems.de/index.php/Tutorial:_Rectangular_Waveguide
%
% Tested with
%  - Matlab 2011a / Octave 3.4.3
%  - openEMS v0.0.26
%
% (C) 2010-2012 Thorsten Liebig <thorsten.liebig@gmx.de>

close all
clear
clc

%% setup the simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
physical_constants;
unit = 1e-3; %drawing unit in mm
numTS = 50000; %max. number of timesteps

% waveguide dimensions
length = 5000;
a = 1000;   %waveguide width
b = 1000;    %waveguide heigth

% frequency range of interest
f_start =  300e6;
f_stop  =  500e6;

%waveguide TE-mode definition
m = 1;
n = 1;

mesh_res = [10 10 10]; %targeted mesh resolution

%% mode functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% by David M. Pozar, Microwave Engineering, third edition, page 113
freq = linspace(f_start,f_stop,201);
k = 2*pi*freq/c0;
kc = sqrt((m*pi/a/unit)^2 + (n*pi/b/unit)^2);
fc = c0*kc/2/pi;          %cut-off frequency
beta = sqrt(k.^2 - kc^2); %waveguide phase-constant
ZL_a = k * Z0 ./ beta;    %analytic waveguide impedance

% mode profile E- and H-field
func_Ex = [num2str( n/b/unit) '*cos(' num2str(m*pi/a) '*x)*sin('  num2str(n*pi/b) '*y)'];
func_Ey = [num2str(-m/a/unit) '*sin(' num2str(m*pi/a) '*x)*cos('  num2str(n*pi/b) '*y)'];

func_Hx = [num2str(m/a/unit) '*sin(' num2str(m*pi/a) '*x)*cos('  num2str(n*pi/b) '*y)'];
func_Hy = [num2str(n/b/unit) '*cos(' num2str(m*pi/a) '*x)*sin('  num2str(n*pi/b) '*y)'];

disp([' Cutoff frequencies for this mode and wavguide is: ' num2str(fc/1e6) ' MHz']);

if (f_start<fc)
    warning('openEMS:example','f_start is smaller than the cutoff-frequency, this may result in a long simulation... ');
end

%% setup FDTD parameter & excitation function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FDTD = InitFDTD(numTS,1e-5);
FDTD = SetGaussExcite(FDTD,0.5*(f_start+f_stop),0.5*(f_stop-f_start));

% boundary conditions
BC = [0 0 0 0 3 3]; %pml in pos. and neg. z-direction
FDTD = SetBoundaryCond(FDTD,BC);

%% setup CSXCAD mesh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CSX = InitCSX();
mesh.x = SmoothMeshLines([0 a], mesh_res(1));
mesh.y = SmoothMeshLines([0 b], mesh_res(2));
mesh.z = SmoothMeshLines([0 length], mesh_res(3));
CSX = DefineRectGrid(CSX, unit,mesh);

%% apply the excitation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% xy-mode profile excitation located directly on top of pml (first 8 z-lines)
CSX = AddExcitation(CSX,'excite',0,[1 1 0]);
weight{1} = func_Ex;
weight{2} = func_Ey;
weight{3} = 0;
CSX = SetExcitationWeight(CSX,'excite',weight);
start=[mesh.x(1)   mesh.y(1)   mesh.z(8) ];
stop =[mesh.x(end) mesh.y(end) mesh.z(8) ];
CSX = AddBox(CSX,'excite',0 ,start,stop);

%% voltage and current definitions using the mode matching probes %%%%%%%%%
%port 1
start = [mesh.x(1)   mesh.y(1)   mesh.z(15)];
stop  = [mesh.x(end) mesh.y(end) mesh.z(15)];
CSX = AddProbe(CSX, 'ut1', 10, 1, [], 'ModeFunction',{func_Ex,func_Ey,0});
CSX = AddBox(CSX,  'ut1',  0 ,start,stop);
CSX = AddProbe(CSX,'it1', 11, 1, [], 'ModeFunction',{func_Hx,func_Hy,0});
CSX = AddBox(CSX,'it1', 0 ,start,stop);

%port 2
start = [mesh.x(1)   mesh.y(1)   mesh.z(end-14)];
stop  = [mesh.x(end) mesh.y(end) mesh.z(end-14)];
CSX = AddProbe(CSX, 'ut2', 10, 1, [], 'ModeFunction',{func_Ex,func_Ey,0});
CSX = AddBox(CSX,  'ut2',  0 ,start,stop);
CSX = AddProbe(CSX,'it2', 11, 1, [], 'ModeFunction',{func_Hx,func_Hy,0});
CSX = AddBox(CSX,'it2', 0 ,start,stop);

port_dist = mesh.z(end-14) - mesh.z(15);
 
%% define dump box... %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CSX = AddDump(CSX,'Et','FileType',1,'SubSampling','4,4,4');
start = [mesh.x(1)   mesh.y(1)   mesh.z(1)];
stop  = [mesh.x(end) mesh.y(end) mesh.z(end)];
CSX = AddBox(CSX,'Et',0 , start,stop);

%% Write openEMS compatoble xml-file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sim_Path = 'tmp';
Sim_CSX = 'rect_wg.xml';

[status, message, messageid] = rmdir(Sim_Path,'s');
[status, message, messageid] = mkdir(Sim_Path);

WriteOpenEMS([Sim_Path '/' Sim_CSX],FDTD,CSX);

RunOpenEMS(Sim_Path, Sim_CSX)

%% postproc %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U = ReadUI({'ut1','ut2'},[Sim_Path '/'],freq);
I = ReadUI({'it1','it2'},[Sim_Path '/'],freq);
Exc = ReadUI('et',Sim_Path,freq);

uf1 = U.FD{1}.val./Exc.FD{1}.val;
uf2 = U.FD{2}.val./Exc.FD{1}.val;
if1 = I.FD{1}.val./Exc.FD{1}.val;
if2 = I.FD{2}.val./Exc.FD{1}.val;

uf1_inc = 0.5 * ( uf1 + if1 .* ZL_a );
if1_inc = 0.5 * ( if1 + uf1 ./ ZL_a );
uf2_inc = 0.5 * ( uf2 + if2 .* ZL_a );
if2_inc = 0.5 * ( if2 + uf2 ./ ZL_a );

uf1_ref = uf1 - uf1_inc;
if1_ref = if1_inc - if1; 
uf2_ref = uf2 - uf2_inc;
if2_ref = if2_inc - if2;

%% plot s-parameter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
s11 = uf1_ref./uf1_inc;
s21 = uf2_inc./uf1_inc;
plot(freq*1e-6,20*log10(abs(s11)),'k-','Linewidth',2);
xlim([freq(1) freq(end)]*1e-6);
grid on;
hold on;
plot(freq*1e-6,20*log10(abs(s21)),'r--','Linewidth',2);
l = legend('S_{11}','S_{21}','Location','Best');
set(l,'FontSize',12);
ylabel('S-Parameter (dB)','FontSize',12);
xlabel('frequency (MHz) \rightarrow','FontSize',12);

%% compare analytic and numerical wave-impedance %%%%%%%%%%%%%%%%%%%%%%%%%%
figure
ZL = uf1./if1;
plot(freq*1e-6,real(ZL),'Linewidth',2);
hold on;
grid on;
plot(freq*1e-6,imag(ZL),'r--','Linewidth',2);
plot(freq*1e-6,ZL_a,'g-.','Linewidth',2);
ylabel('ZL (\Omega)','FontSize',12);
xlabel('frequency (MHz) \rightarrow','FontSize',12);
xlim([freq(1) freq(end)]*1e-6);
l = legend('\Re(Z_L)','\Im(Z_L)','Z_L analytic','Location','Best');
set(l,'FontSize',12);

%% Plot the field dumps %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
dump_file = [Sim_Path '/Et.h5'];
PlotArgs.slice = {a/2*unit b/2*unit 0};
PlotArgs.pauseTime=0.01;
PlotArgs.component=0;
PlotArgs.Limit = 'auto';
PlotHDF5FieldData(dump_file, PlotArgs)
