%
% Tutorials / Circ_Waveguide
%
% Describtion at:
% http://openems.de/index.php/Tutorial:_Circular_Waveguide
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

% waveguide dimensions
length = 2000;
rad = 350;     %waveguide radius in mm

% frequency range of interest
f_start =  300e6;
f_stop  =  500e6;

mesh_res = [10 2*pi/49.999 10]; %targeted mesh resolution

%% mode functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% by David M. Pozar, Microwave Engineering, third edition
freq = linspace(f_start,f_stop,201);
p11 = 1.841;
kc = p11 / rad /unit;
k = 2*pi*freq/C0;
fc = C0*kc/2/pi;
beta = sqrt(k.^2 - kc^2);
n_eff = (beta/k);
ZL_a = k * Z0 ./ beta;    %analytic waveguide impedance

% TE_11 mode profile E- and H-field
kc = kc*unit; %functions must be defined in drawing units
func_Er = [ num2str(-1/kc^2,15) '/rho*cos(a)*j1('  num2str(kc,15) '*rho)'];
func_Ea = [ num2str(1/kc,15) '*sin(a)*0.5*(j0('  num2str(kc,15) '*rho)-jn(2,'  num2str(kc,15) '*rho))'];

func_Ha = [ num2str(-1/kc^2,15) '/rho*cos(a)*j1('  num2str(kc,15) '*rho)'];
func_Hr = [ '-1*' num2str(1/kc,15) '*sin(a)*0.5*(j0('  num2str(kc,15) '*rho)-jn(2,'  num2str(kc,15) '*rho))'];

disp([' Cutoff frequencies for this mode and wavguide is: ' num2str(fc/1e6) ' MHz']);

if (f_start<fc)
    warning('openEMS:example','f_start is smaller than the cutoff-frequency, this may result in a long simulation... ');
end

%% setup FDTD parameter & excitation function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FDTD = InitFDTD(1e6,1e-5,'CoordSystem',1);
FDTD = SetGaussExcite(FDTD,0.5*(f_start+f_stop),0.5*(f_stop-f_start));

% boundary conditions
BC = [0 0 0 0 3 3]; %pml in pos. and neg. z-direction
FDTD = SetBoundaryCond(FDTD,BC);

%% setup CSXCAD mesh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CSX = InitCSX('CoordSystem',1); % init a cylindrical mesh
mesh.r = SmoothMeshLines([0 rad], mesh_res(1)); %mesh in radial direction
mesh.a = SmoothMeshLines([0 2*pi], mesh_res(2)); % mesh in aziumthal dir.
mesh.z = SmoothMeshLines([0 length], mesh_res(3));
CSX = DefineRectGrid(CSX, unit,mesh);

%% apply the excitation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ra-mode profile excitation located directly on top of pml (first 8 z-lines)
CSX = AddExcitation(CSX,'excite',0,[1 1 0]);
weight{1} = func_Er;
weight{2} = func_Ea;
weight{3} = 0;
CSX = SetExcitationWeight(CSX,'excite',weight);
start=[mesh.r(1)   mesh.a(1)   mesh.z(8) ];
stop =[mesh.r(end) mesh.a(end) mesh.z(8) ];
CSX = AddBox(CSX,'excite',0 ,start,stop);

%% voltage and current definitions using the mode matching probes %%%%%%%%%
%port 1
start = [mesh.r(1)   mesh.a(1)   mesh.z(15)];
stop  = [mesh.r(end) mesh.a(end) mesh.z(15)];
CSX = AddProbe(CSX, 'ut1', 10, 'ModeFunction',{func_Er,func_Ea,0});
CSX = AddBox(CSX,  'ut1',  0 ,start,stop);
CSX = AddProbe(CSX,'it1', 11, 'ModeFunction',{func_Hr,func_Ha,0});
CSX = AddBox(CSX,'it1', 0 ,start,stop);

%port 2
start(3) = mesh.z(end-14);
stop(3) = mesh.z(end-14);
CSX = AddProbe(CSX, 'ut2', 10, 'ModeFunction',{func_Er,func_Ea,0});
CSX = AddBox(CSX,  'ut2',  0 ,start,stop);
CSX = AddProbe(CSX,'it2', 11, 'ModeFunction',{func_Hr,func_Ha,0});
CSX = AddBox(CSX,'it2', 0 ,start,stop);

port_dist = mesh.z(end-14) - mesh.z(15);
 
%% define dump box... %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CSX = AddDump(CSX,'Et','FileType',1,'SubSampling','4,4,4');
start = [mesh.r(1)   mesh.a(1)   mesh.z(1)];
stop  = [mesh.r(end) mesh.a(end) mesh.z(end)];
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
