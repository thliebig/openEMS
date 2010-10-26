close all
clear
clc

%% setup the simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
length = 2000;
unit = 1e-3;
a = 1000;
width = a;
b = 500;
height = b;
mesh_res = [10 10 10];

%define mode
m = 1;
n = 0;

EPS0 = 8.85418781762e-12;
MUE0 = 1.256637062e-6;
C0 = 1/sqrt(EPS0*MUE0);
Z0 = sqrt(MUE0/EPS0);

f0 = 1e9;
freq = linspace(f0-f0/3,f0+f0/3,201);

k = 2*pi*freq/C0;
kc = sqrt((m*pi/a/unit)^2 + (n*pi/b/unit)^2);
fc = C0*kc/2/pi;
beta = sqrt(k.^2 - kc^2);

ZL_a = k * Z0 ./ beta; 

%% mode functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
func_Ex = [num2str( n/b/unit) '*cos(' num2str(m*pi/a) '*x)*sin('  num2str(n*pi/b) '*y)'];
func_Ey = [num2str(-m/a/unit) '*sin(' num2str(m*pi/a) '*x)*cos('  num2str(n*pi/b) '*y)'];

func_Hx = [num2str(m/a/unit) '*sin(' num2str(m*pi/a) '*x)*cos('  num2str(n*pi/b) '*y)'];
func_Hy = [num2str(n/b/unit) '*cos(' num2str(m*pi/a) '*x)*sin('  num2str(n*pi/b) '*y)'];

%% define and openEMS options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
openEMS_opts = '';
% openEMS_opts = [openEMS_opts ' --disable-dumps'];
% openEMS_opts = [openEMS_opts ' --debug-material'];
openEMS_opts = [openEMS_opts ' --engine=fastest'];

Settings = [];
Settings.LogFile = 'openEMS.log';

Sim_Path = 'tmp';
Sim_CSX = 'rect_wg.xml';

if (exist(Sim_Path,'dir'))
    rmdir(Sim_Path,'s');
end
mkdir(Sim_Path);

%% setup FDTD parameter & excitation function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FDTD = InitFDTD(50000,1e-5,'OverSampling',6);
FDTD = SetGaussExcite(FDTD,f0,f0/3);
BC = [0 0 0 0 0 3];
FDTD = SetBoundaryCond(FDTD,BC);

%% setup CSXCAD geometry & mesh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CSX = InitCSX();
mesh.x = 0 : mesh_res(1) : width;
mesh.y = 0 : mesh_res(2) : height;
mesh.z = 0 : mesh_res(3) : length;
CSX = DefineRectGrid(CSX, unit,mesh);

%% apply the excitation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
start=[0     0      mesh.z(1) ];
stop =[width height mesh.z(1) ];
CSX = AddExcitation(CSX,'excite',0,[1 1 0]);
weight{1} = func_Ex;
weight{2} = func_Ey;
weight{3} = 0;
CSX = SetExcitationWeight(CSX,'excite',weight);
CSX = AddBox(CSX,'excite',0 ,start,stop);

%% voltage and current definitions using the mode matching probes %%%%%%%%%
start = [mesh.x(1)   mesh.y(1)   mesh.z(15)];
stop  = [mesh.x(end) mesh.y(end) mesh.z(15)];
CSX = AddProbe(CSX, 'ut1', 10, 1, [], 'ModeFunction',{func_Ex,func_Ey,0});
CSX = AddBox(CSX,  'ut1',  0 ,start,stop);
CSX = AddProbe(CSX,'it1', 11, 1, [], 'ModeFunction',{func_Hx,func_Hy,0});
CSX = AddBox(CSX,'it1', 0 ,start,stop);

start = [mesh.x(1)   mesh.y(1)   mesh.z(end-15)];
stop  = [mesh.x(end) mesh.y(end) mesh.z(end-15)];
CSX = AddProbe(CSX, 'ut2', 10, 1, [], 'ModeFunction',{func_Ex,func_Ey,0});
CSX = AddBox(CSX,  'ut2',  0 ,start,stop);
CSX = AddProbe(CSX,'it2', 11, 1, [], 'ModeFunction',{func_Hx,func_Hy,0});
CSX = AddBox(CSX,'it2', 0 ,start,stop);
 
%% define dump boxes... %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CSX = AddDump(CSX,'Et','FileType',1,'SubSampling','4,4,2');
start = [mesh.x(1) , height/2 , mesh.z(1)];
stop = [mesh.x(end) ,  height/2 , mesh.z(end)];
CSX = AddBox(CSX,'Et',0 , start,stop);

CSX = AddDump(CSX,'Ht','DumpType',1,'FileType',1,'SubSampling','4,4,2');
CSX = AddBox(CSX,'Ht',0,start,stop);

%% Write openEMS compatoble xml-file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WriteOpenEMS([Sim_Path '/' Sim_CSX],FDTD,CSX);

RunOpenEMS(Sim_Path, Sim_CSX, openEMS_opts, Settings)


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

%% plot s-parameter
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

%% compare analytic and numerical wave-impedance
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

