%
% EXAMPLE / waveguide / circular waveguide
%
% This example demonstrates how to:
%  - setup a circular waveguide
%  - use analytic functions for waveguide excitations and voltage/current
%  calculations
% 
%
% Tested with
%  - Matlab 2009b
%  - openEMS v0.0.17
%
% (C) 2010 Thorsten Liebig <thorsten.liebig@uni-due.de>

close all
clear
clc

%% switches & options...
postprocessing_only = 0;
use_pml = 0;            % use pml boundaries instead of mur
openEMS_opts = '';
% openEMS_opts = [openEMS_opts ' --disable-dumps'];

%% setup the simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numTS = 1e5;            %number of timesteps
length = 1000;          %length of the waveguide
unit = 1e-3;            %drawing unit used
rad  = 300;             %radius of the circular waveguide
mesh_res = [10 10 15];  %desired mesh resolution

%excitation
f0 = 350e6;             %center frequency
f0_BW = 25e6;           %bandwidth: 10dB cut-off frequency

physical_constants

%% TE11 mode definitions (Pozar 3rd edition)
p11 = 1.841;
kc = p11 / rad /unit;
k = 2*pi*f0/C0;
fc = C0*kc/2/pi;
beta = sqrt(k^2 - kc^2);
n_eff = (beta/k);

kc = kc*unit; %functions must be defined in drawing units
func_Er = [ num2str(-1/kc^2) '/rho*cos(a)*j1('  num2str(kc) '*rho)'];
func_Ea = [ num2str(1/kc) '*sin(a)*0.5*(j0('  num2str(kc) '*rho)-jn(2,'  num2str(kc) '*rho))'];
func_Ex = ['(' func_Er '*cos(a) - ' func_Ea '*sin(a) )*(rho<' num2str(rad) ')'];
func_Ey = ['(' func_Er '*sin(a) + ' func_Ea '*cos(a) )*(rho<' num2str(rad) ')'];

func_Ha = [ num2str(-1/kc^2,'%14.13f') '/rho*cos(a)*j1('  num2str(kc,'%14.13f') '*rho)'];
func_Hr = [ '-1*' num2str(1/kc,'%14.13f') '*sin(a)*0.5*(j0('  num2str(kc,'%14.13f') '*rho)-jn(2,'  num2str(kc,'%14.13f') '*rho))'];
func_Hx = ['(' func_Hr '*cos(a) - ' func_Ha '*sin(a) )*(rho<' num2str(rad) ')'];
func_Hy = ['(' func_Hr '*sin(a) + ' func_Ha '*cos(a) )*(rho<' num2str(rad) ')'];

%% define files and path %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sim_Path = 'tmp';
Sim_CSX = 'Circ_WG.xml';

if (postprocessing_only==0)
    [status, message, messageid] = rmdir(Sim_Path,'s');
    [status, message, messageid] = mkdir(Sim_Path);
end

%% setup FDTD parameter & excitation function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FDTD = InitFDTD(numTS,1e-6,'OverSampling',5);
FDTD = SetGaussExcite(FDTD,f0,f0_BW);
BC = {'PEC','PEC','PEC','PEC','PEC','MUR'};
if (use_pml>0)
    BC = {'PEC','PEC','PEC','PEC','PEC','PML_8'};
end
FDTD = SetBoundaryCond(FDTD,BC,'MUR_PhaseVelocity',C0 / n_eff);

%% setup CSXCAD geometry & mesh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CSX = InitCSX();
mesh.x = -mesh_res(1)/2-rad:mesh_res(1):rad+mesh_res(1)/2;
mesh.y = -mesh_res(2)/2-rad:mesh_res(2):rad+mesh_res(2)/2;
mesh.z = 0 : mesh_res(3) : length;
CSX = DefineRectGrid(CSX, 1e-3,mesh);

start = [0,0,0];
stop = [0,0,length];

%%% fill everything with copper, priority 0
CSX = AddMetal(CSX,'copper');
% CSX = SetMaterialProperty(CSX,'copper','Kappa',56e6);
CSX = AddBox(CSX,'copper',0,[mesh.x(1) mesh.y(1) mesh.z(1)],[mesh.x(end) mesh.y(end) mesh.z(end)]);

%%% cut out an air cylinder as circular waveguide... priority 5
CSX = AddMaterial(CSX,'air');
CSX = SetMaterialProperty(CSX,'air','Epsilon',1);
CSX = AddCylinder(CSX,'air', 5 ,start,stop,rad);

CSX = AddExcitation(CSX,'excite',0,[1 1 0]);
weight{1} = func_Ex;
weight{2} = func_Ey;
weight{3} = 0;
CSX = SetExcitationWeight(CSX, 'excite', weight );
CSX = AddCylinder(CSX,'excite', 5 ,[0 0 -0.1],[0 0 0.1],rad);
 
%% define dump boxes... %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CSX = AddDump(CSX,'Et_','SubSampling','2,2,2','FileType',0,'DumpMode',2);
start = [mesh.x(1) , 0 , mesh.z(1)];
stop = [mesh.x(end), 0 , mesh.z(end)];
CSX = AddBox(CSX,'Et_',0 , start,stop);

CSX = AddDump(CSX,'Ht_','SubSampling','2,2,2','DumpType',1,'FileType',0,'DumpMode',2);
CSX = AddBox(CSX,'Ht_',0,start,stop);

%% define voltage calc boxes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%voltage calc
start = [mesh.x(1)   mesh.y(1)   mesh.z(10)];
stop  = [mesh.x(end) mesh.y(end) mesh.z(10)];
CSX = AddProbe(CSX, 'ut1', 10, 1, [], 'ModeFunction',{func_Ex,func_Ey,0});
CSX = AddBox(CSX,  'ut1',  0 ,start,stop);
CSX = AddProbe(CSX,'it1', 11, 1, [], 'ModeFunction',{func_Hx,func_Hy,0});
CSX = AddBox(CSX,'it1', 0 ,start,stop);
    
start = [mesh.x(1)   mesh.y(1)   mesh.z(end-10)];
stop  = [mesh.x(end) mesh.y(end) mesh.z(end-10)];
CSX = AddProbe(CSX, 'ut2', 10, 1, [], 'ModeFunction',{func_Ex,func_Ey,0});
CSX = AddBox(CSX,  'ut2',  0 ,start,stop);
CSX = AddProbe(CSX,'it2', 11, 1, [], 'ModeFunction',{func_Hx,func_Hy,0});
CSX = AddBox(CSX,'it2', 0 ,start,stop);

port_dist = mesh.z(end-10) - mesh.z(10);

%% Write openEMS
if (postprocessing_only==0)
    WriteOpenEMS([Sim_Path '/' Sim_CSX],FDTD,CSX);

    RunOpenEMS(Sim_Path, Sim_CSX, openEMS_opts);
end

%% do the plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
freq = linspace(f0-f0_BW,f0+f0_BW,201);
U = ReadUI({'ut1','ut2'},[Sim_Path '/'],freq);
I = ReadUI({'it1','it2'},[Sim_Path '/'],freq);
Exc = ReadUI('et',Sim_Path,freq);

k = 2*pi*freq/C0;
kc = p11 / rad /unit;
beta = sqrt(k.^2 - kc^2);

ZL_a = Z0*k./beta ;

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

% plot s-parameter
figure
s11 = uf1_ref./uf1_inc;
s21 = uf2_inc./uf1_inc;
plot(freq,20*log10(abs(s11)),'Linewidth',2);
xlim([freq(1) freq(end)]);
xlabel('frequency (Hz)')
ylabel('s-para (dB)');
% ylim([-40 5]);
grid on;
hold on;
plot(freq,20*log10(abs(s21)),'r','Linewidth',2);
legend('s11','s21','Location','SouthEast');

% plot line-impedance comparison
figure()
ZL = uf1./if1;
plot(freq,real(ZL),'Linewidth',2);
xlim([freq(1) freq(end)]);
xlabel('frequency (Hz)')
ylabel('line-impedance (\Omega)');
grid on;
hold on;
plot(freq,imag(ZL),'r--','Linewidth',2);
plot(freq,ZL_a,'g-.','Linewidth',2);
legend('\Re\{ZL\}','\Im\{ZL\}','ZL-analytic','Location','Best');

% beta compare
figure()
da = angle(uf1_inc)-angle(uf2_inc);
da = mod(da,2*pi);
beta_12 = (da)/port_dist/unit;
plot(freq,beta_12,'Linewidth',2);
xlim([freq(1) freq(end)]);
xlabel('frequency (Hz)');
ylabel('\beta (m^{-1})');
grid on;
hold on;
plot(freq,beta,'g--','Linewidth',2);
legend('\beta-FDTD','\beta-analytic','Location','Best');

%% visualize electric & magnetic fields
disp('you will find vtk dump files in the simulation folder (tmp/)')
disp('use paraview to visulaize them');
