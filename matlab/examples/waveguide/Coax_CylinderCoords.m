%
% EXAMPLE / waveguide / coaxial cable using cylindrical coordinates
%
% This example demonstrates how to:
%  - use cylindrical coordinates
%  - setup a coaxial waveguide
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
coax_rad_i  = 100;      %inner radius
coax_rad_a  = 230;      %outer radius
mesh_res = [10 nan 10]; %desired mesh resolution
N_alpha = 71;           %mesh lines in azimuth direction

physical_constants;

%excitation
f0 = 0.5e9;
epsR = 1;

%% create sim path %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sim_Path = 'tmp';
Sim_CSX = 'coax.xml';

if (postprocessing_only==0)
    [status, message, messageid] = rmdir(Sim_Path,'s');
    [status, message, messageid] = mkdir(Sim_Path);
end

%% setup FDTD parameter & excitation function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FDTD = InitCylindricalFDTD(numTS,1e-5);
FDTD = SetGaussExcite(FDTD,f0,f0);
BC = {'PEC','PEC','PEC','PEC','PEC','MUR'};
if (use_pml>0)
    BC = {'PEC','PEC','PEC','PEC','PEC','PML_8'};
end
FDTD = SetBoundaryCond(FDTD,BC);

%% setup CSXCAD geometry & mesh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CSX = InitCSX('CoordSystem',1);
mesh.x = coax_rad_i : mesh_res(1) : coax_rad_a;
mesh.y = linspace(0,2*pi,N_alpha);
mesh.z = 0 : mesh_res(3) : length;
CSX = DefineRectGrid(CSX, unit, mesh);

%% material
CSX = AddMaterial(CSX,'fill');
CSX = SetMaterialProperty(CSX,'fill','Epsilon',epsR);
start = [mesh.x(1) mesh.y(1) 0];
stop = [mesh.x(end) mesh.y(end) length];
CSX = AddBox(CSX,'fill',0 ,start,stop);

%% apply the excitation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CSX = AddExcitation(CSX,'excite',0,[1 0 0]);
weight{1} = '1/rho';
weight{2} = 0;
weight{3} = 0;
CSX = SetExcitationWeight(CSX, 'excite', weight );
start = [coax_rad_i mesh.y(1)   0];
stop  = [coax_rad_a mesh.y(end) 0];
CSX = AddBox(CSX,'excite',0 ,start,stop);
 
%% define dump boxes... %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CSX = AddDump(CSX,'Et_','DumpMode',0);
start = [mesh.x(1) , 0 , mesh.z(1)];
stop = [mesh.x(end) , 0 , mesh.z(end)];
CSX = AddBox(CSX,'Et_',0 , start,stop);

CSX = AddDump(CSX,'Ht_','DumpType',1,'DumpMode',0);
CSX = AddBox(CSX,'Ht_',0,start,stop);

%% define voltage calc boxes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%voltage calc
CSX = AddProbe(CSX,'ut1',0);
start = [ coax_rad_i 0 mesh.z(10) ];
stop  = [ coax_rad_a 0 mesh.z(10) ];
CSX = AddBox(CSX,'ut1', 0 ,start,stop);
CSX = AddProbe(CSX,'ut2',0);
start = [ coax_rad_i 0 mesh.z(end-10)];
stop  = [ coax_rad_a 0 mesh.z(end-10)];
CSX = AddBox(CSX,'ut2', 0 ,start,stop);

%current calc, for each position there are two currents, which will get
%averaged to match the voltage position in between (!Yee grid!)
CSX = AddProbe(CSX,'it1a',1);
mid = 0.5*(coax_rad_i+coax_rad_a);
start = [ 0   mesh.z(1)   mesh.z(9) ];
stop  = [ mid mesh.z(end) mesh.z(9) ];
CSX = AddBox(CSX,'it1a', 0 ,start,stop);
CSX = AddProbe(CSX,'it1b',1);
start = [ 0   mesh.z(1)   mesh.z(10) ];
stop  = [ mid mesh.z(end) mesh.z(10) ];
CSX = AddBox(CSX,'it1b', 0 ,start,stop);

CSX = AddProbe(CSX,'it2a',1);
start = [ 0   mesh.z(1)   mesh.z(end-11) ];
stop  = [ mid mesh.z(end) mesh.z(end-11) ];
CSX = AddBox(CSX,'it2a', 0 ,start,stop);
CSX = AddProbe(CSX,'it2b',1);
start = [ 0   mesh.z(1)   mesh.z(end-10) ];
stop  = [ mid mesh.z(end) mesh.z(end-10) ];
CSX = AddBox(CSX,'it2b', 0 ,start,stop);

%% Write openEMS
if (postprocessing_only==0)
    WriteOpenEMS([Sim_Path '/' Sim_CSX],FDTD,CSX);
    RunOpenEMS(Sim_Path, Sim_CSX, openEMS_opts);
end

%%
freq = linspace(0,2*f0,201);
U = ReadUI({'ut1','ut2'},[Sim_Path '/'],freq);
I = ReadUI({'it1a','it1b','it2a','it2b'},[Sim_Path '/'],freq);
Exc = ReadUI('et',Sim_Path,freq);

%% plot voltages
figure
plot(U.TD{1}.t, U.TD{1}.val,'Linewidth',2);
hold on;
grid on;
plot(U.TD{2}.t, U.TD{2}.val,'r--','Linewidth',2);
xlabel('time (s)')
ylabel('voltage (V)')
legend('u_1(t)','u_2(t)')

%% calculate incoming and reflected voltages & currents
ZL_a = ones(size(freq))*Z0/2/pi/sqrt(epsR)*log(coax_rad_a/coax_rad_i); %analytic line-impedance of a coax

uf1 = U.FD{1}.val./Exc.FD{1}.val;
uf2 = U.FD{2}.val./Exc.FD{1}.val;
if1 = 0.5*(I.FD{1}.val+I.FD{2}.val)./Exc.FD{1}.val;
if2 = 0.5*(I.FD{3}.val+I.FD{4}.val)./Exc.FD{1}.val;

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
