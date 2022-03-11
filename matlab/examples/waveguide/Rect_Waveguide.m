%
% EXAMPLE / waveguide / Rect_Waveguide
%
% This example demonstrates:
%  - waveguide mode excitation
%  - waveguide mode matching
%  - pml absorbing boundaries
% 
%
% Tested with
%  - Matlab 2009b
%  - openEMS v0.0.17
%
% (C) 2010 Thorsten Liebig <thorsten.liebig@gmx.de>

close all
clear
clc

%% switches
postproc_only = 0;

%% setup the simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
physical_constants;
unit = 1e-3; %drawing unit in mm
numTS = 50000; %max. number of timesteps

% waveguide dimensions
length = 1000;
a = 1000;   %waveguide width
b = 600;    %waveguide height

%waveguide TE-mode definition
m = 1;
n = 0;

mesh_res = [10 10 10];

%% setup FDTD parameters & excitation function %%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_start =  175e6;
f_stop  =  500e6;

% dump special frequencies to vtk, use paraview (www.paraview.org) to
% animate this dumps over phase
vtk_dump_freq = [200e6 300e6 500e6];

freq = linspace(f_start,f_stop,201);

k = 2*pi*freq/c0;
kc = sqrt((m*pi/a/unit)^2 + (n*pi/b/unit)^2);
fc = c0*kc/2/pi;          %cut-off frequency
beta = sqrt(k.^2 - kc^2); %waveguide phase-constant
ZL_a = k * Z0 ./ beta;    %analytic waveguide impedance

disp([' Cutoff frequencies for this mode and wavguide is: ' num2str(fc/1e6) ' MHz']);

if (f_start<fc)
    warning('openEMS:example','f_start is smaller than the cutoff-frequency, this may result in a long simulation... ');
end

%% mode functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% by David M. Pozar, Microwave Engineering, third edition, page 113
func_Ex = [num2str( n/b/unit) '*cos(' num2str(m*pi/a) '*x)*sin('  num2str(n*pi/b) '*y)'];
func_Ey = [num2str(-m/a/unit) '*sin(' num2str(m*pi/a) '*x)*cos('  num2str(n*pi/b) '*y)'];

func_Hx = [num2str(m/a/unit) '*sin(' num2str(m*pi/a) '*x)*cos('  num2str(n*pi/b) '*y)'];
func_Hy = [num2str(n/b/unit) '*cos(' num2str(m*pi/a) '*x)*sin('  num2str(n*pi/b) '*y)'];

%% define and openEMS options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
openEMS_opts = '';
% openEMS_opts = [openEMS_opts ' --disable-dumps'];
% openEMS_opts = [openEMS_opts ' --debug-material'];
% openEMS_opts = [openEMS_opts ' --engine=basic'];

Settings = [];
Settings.LogFile = 'openEMS.log';

Sim_Path = 'tmp';
Sim_CSX = 'rect_wg.xml';

if (postproc_only==0)
    [status, message, messageid] = rmdir(Sim_Path,'s');
    [status, message, messageid] = mkdir(Sim_Path);
end

%% setup FDTD parameter & excitation function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FDTD = InitFDTD(numTS,1e-5,'OverSampling',6);
FDTD = SetGaussExcite(FDTD,0.5*(f_start+f_stop),0.5*(f_stop-f_start));
BC = [0 0 0 0 0 3];
FDTD = SetBoundaryCond(FDTD,BC);

%% setup CSXCAD geometry & mesh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CSX = InitCSX();
mesh.x = SmoothMeshLines([0 a], mesh_res(1));
mesh.y = SmoothMeshLines([0 b], mesh_res(2));
mesh.z = SmoothMeshLines([0 length], mesh_res(3));
CSX = DefineRectGrid(CSX, unit,mesh);

%% apply the excitation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
start=[mesh.x(1)   mesh.y(1)   mesh.z(1) ];
stop =[mesh.x(end) mesh.y(end) mesh.z(1) ];
CSX = AddExcitation(CSX,'excite',0,[1 1 0]);
weight{1} = func_Ex;
weight{2} = func_Ey;
weight{3} = 0;
CSX = SetExcitationWeight(CSX,'excite',weight);
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
start = [mesh.x(1)   mesh.y(1)   mesh.z(end-15)];
stop  = [mesh.x(end) mesh.y(end) mesh.z(end-15)];
CSX = AddProbe(CSX, 'ut2', 10, 1, [], 'ModeFunction',{func_Ex,func_Ey,0});
CSX = AddBox(CSX,  'ut2',  0 ,start,stop);
CSX = AddProbe(CSX,'it2', 11, 1, [], 'ModeFunction',{func_Hx,func_Hy,0});
CSX = AddBox(CSX,'it2', 0 ,start,stop);

port_dist = mesh.z(end-15) - mesh.z(15);
 
%% define dump boxes... %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CSX = AddDump(CSX,'Et','FileType',1,'SubSampling','4,4,2');
start = [mesh.x(1)   mesh.y(1)   mesh.z(1)];
stop  = [mesh.x(end) mesh.y(end) mesh.z(end)];
CSX = AddBox(CSX,'Et',0 , start,stop);

CSX = AddDump(CSX,'Ht','DumpType',1,'FileType',1,'SubSampling','4,4,2');
CSX = AddBox(CSX,'Ht',0,start,stop);

%% Write openEMS compatoble xml-file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (postproc_only==0)
    WriteOpenEMS([Sim_Path '/' Sim_CSX],FDTD,CSX);

    RunOpenEMS(Sim_Path, Sim_CSX, openEMS_opts, Settings)
end

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
if1_ref = if1 - if1_inc;
uf2_ref = uf2 - uf2_inc;
if2_ref = if2 - if2_inc;

%% plot s-parameter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
s11 = uf1_ref./uf1_inc;
s21 = uf2_inc./uf1_inc;
plot(freq,20*log10(abs(s11)),'Linewidth',2);
xlim([freq(1) freq(end)]);
% ylim([-40 5]);
grid on;
hold on;
plot(freq,20*log10(abs(s21)),'r','Linewidth',2);
legend('s11','s21','Location','SouthEast');
ylabel('s-para (dB)');
xlabel('freq (Hz)');

%% compare analytic and numerical wave-impedance %%%%%%%%%%%%%%%%%%%%%%%%%%
ZL = uf1./if1;
figure()
plot(freq,real(ZL),'Linewidth',2);
hold on;
grid on;
plot(freq,imag(ZL),'r--','Linewidth',2);
plot(freq,ZL_a,'g-.','Linewidth',2);
ylabel('ZL (\Omega)');
xlabel('freq (Hz)');
xlim([freq(1) freq(end)]);
legend('\Re(Z_L)','\Im(Z_L)','Z_L analytic','Location','Best');

%% beta compare
figure()
da = unwrap(angle(uf1_inc./uf2_inc)) ;
% da = mod(da,2*pi);
beta_12 = (da)/port_dist/unit;
plot(freq,beta_12,'Linewidth',2);
xlim([freq(1) freq(end)]);
xlabel('frequency (Hz)');
ylabel('\beta (m^{-1})');
grid on;
hold on;
plot(freq,beta,'g--','Linewidth',2);
legend('\beta-FDTD','\beta-analytic','Location','Best');

%% Plot the field dumps %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dump_file = [Sim_Path '/Et.h5'];
figure()
PlotArgs.slice = {a/2*unit b/2*unit 0};
PlotArgs.pauseTime=0.01;
PlotArgs.component=0;
PlotArgs.Limit = 'auto';
PlotHDF5FieldData(dump_file, PlotArgs)

%% dump frequency to vtk %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cleanup and create dump folder
vtk_path = [Sim_Path '/vtk'];
[status, message, messageid] = rmdir(vtk_path,'s');
[status, message, messageid] = mkdir(vtk_path);

disp('Dumping to vtk files... this may take a minute...')
% define interpolation mesh
mesh_interp{1}=mesh.x * unit;
mesh_interp{2}=b/2 * unit;
mesh_interp{3}=mesh.z * unit;
[field mesh_FD] = ReadHDF5Dump(dump_file,'Interpolation',mesh_interp,'Frequency',vtk_dump_freq);

% dump animated phase to vtk
for n=1:numel(vtk_dump_freq)   
    phase = linspace(0,360,21);
    phase = phase(1:end-1);
    for ph = phase
        filename = [vtk_path '/E_xz_f=' num2str(vtk_dump_freq(n)) '_p' num2str(ph) '.vtk'];
        Dump2VTK(filename,real(field.FD.values{n}.*exp(1j*ph/180*pi)),mesh_FD,'E-Field');
    end
    
    filename = [vtk_path '/E_xz_f=' num2str(vtk_dump_freq(n)) '_mag.vtk'];
    Dump2VTK(filename,abs(field.FD.values{n}),mesh_FD,'E-Field');
end

disp('done... you can open and visualize the vtk-files using Paraview (www.paraview.org)!')
