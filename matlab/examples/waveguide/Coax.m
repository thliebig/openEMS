%
% EXAMPLE / waveguide / coaxial cable
%
% This example demonstrates how to:
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
numTS = 5000;           %number of timesteps
length = 1000;          %length of the waveguide
unit = 1e-3;            %drawing unit used
coax_rad_i  = 100;      %inner radius
coax_rad_ai = 230;      %inner radius of outer cladding
coax_rad_aa = 240;      %outer radius of outer cladding
mesh_res = [5 5 5];     %desired mesh resolution

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
FDTD = InitFDTD(numTS,1e-5);
FDTD = SetGaussExcite(FDTD,f0,f0);
BC = {'PEC','PEC','PEC','PEC','MUR','MUR'};
if (use_pml>0)
    BC = {'PEC','PEC','PEC','PEC','PML_8','PML_8'};
end
FDTD = SetBoundaryCond(FDTD,BC);

%% setup CSXCAD geometry & mesh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CSX = InitCSX();
mesh.x = -coax_rad_aa : mesh_res(1) : coax_rad_aa;
mesh.y = mesh.x;
mesh.z = SmoothMeshLines([0 length], mesh_res(3));
CSX = DefineRectGrid(CSX, unit, mesh);

%%% coax
CSX = AddMetal(CSX,'copper');
start = [0,0,0];
stop  = [0,0,length/2];
[CSX,port{1}] = AddCoaxialPort( CSX, 10, 1, 'copper', '', start, stop, 'z', coax_rad_i, coax_rad_ai, coax_rad_aa, 'ExciteAmp', 1,'FeedShift', 10*mesh_res(1) );

start = [0,0,length/2];
stop  = [0,0,length];
[CSX,port{2}] = AddCoaxialPort( CSX, 10, 2, 'copper', '', start, stop, 'z', coax_rad_i, coax_rad_ai, coax_rad_aa );
 
%% define dump boxes... %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CSX = AddDump(CSX,'Et_','DumpMode',2);
start = [mesh.x(1) , 0 , mesh.z(1)];
stop = [mesh.x(end) , 0 , mesh.z(end)];
CSX = AddBox(CSX,'Et_',0 , start,stop);

CSX = AddDump(CSX,'Ht_','DumpType',1,'DumpMode',2);
CSX = AddBox(CSX,'Ht_',0,start,stop);

%% Write openEMS
if (postprocessing_only==0)
    WriteOpenEMS([Sim_Path '/' Sim_CSX],FDTD,CSX);
    CSXGeomPlot([Sim_Path '/' Sim_CSX]);
    RunOpenEMS(Sim_Path, Sim_CSX, openEMS_opts);
end

%%
freq = linspace(0,2*f0,201);
port = calcPort(port, Sim_Path, freq);

%% plot s-parameter
figure
s11 = port{1}.uf.ref./port{1}.uf.inc;
s21 = port{2}.uf.inc./port{1}.uf.inc;
plot(freq,20*log10(abs(s11)),'Linewidth',2);
hold on
grid on
plot(freq,20*log10(abs(s21)),'r--','Linewidth',2);
xlim([freq(1) freq(end)]);
xlabel('frequency (Hz)')
ylabel('s-para (dB)');

%% plot line-impedance comparison
figure()
ZL_a = ones(size(freq))*Z0/2/pi/sqrt(epsR)*log(coax_rad_ai/coax_rad_i); %analytic line-impedance of a coax
ZL = port{2}.uf.tot./port{2}.if.tot;
plot(freq,real(port{1}.ZL),'Linewidth',2);
xlim([freq(1) freq(end)]);
xlabel('frequency (Hz)')
ylabel('line-impedance (\Omega)');
grid on;
hold on;
plot(freq,imag(port{1}.ZL),'r--','Linewidth',2);
plot(freq,ZL_a,'g-.','Linewidth',2);
legend('\Re\{ZL\}','\Im\{ZL\}','ZL-analytic','Location','Best');
