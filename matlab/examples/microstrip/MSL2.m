%
% microstrip line example
%
% this example shows how to use a MSL port
%
% The MSL is excited at the center of the computational volume. The
% boundary at xmin is an absorbing boundary (Mur) and at xmax an electric
% wall. The reflection coefficient at this wall is S11 = -1.
%


close all
clear
clc

physical_constants


postprocessing_only = 0;


%% setup the simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawingunits = 1e-6; % specify everything in um
MSL_length   = 10000;
MSL_width    = 1000;
substrate_thickness = 254;

mesh_res      = [200 0 0];
max_timesteps = 20000;
min_decrement = 1e-6;
f_max         = 8e9;

%% setup FDTD parameters & excitation function %%%%%%%%%%%%%%%%%%%%%%%%%%%%
FDTD = InitFDTD( max_timesteps, min_decrement, 'OverSampling', 10 );
FDTD = SetGaussExcite( FDTD, f_max/2, f_max/2 );
BC   = [2 0 0 0 0 1];
FDTD = SetBoundaryCond( FDTD, BC );

%% setup CSXCAD geometry & mesh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CSX = InitCSX();
mesh.x = -MSL_length : mesh_res(1) : MSL_length;
mesh.y = linspace(-MSL_width/2,MSL_width/2,10); % discretize the width of the MSL with 10 cells
temp1 = linspace(-4*MSL_width,mesh.y(1),20);
temp2 = linspace(mesh.y(end),4*MSL_width,20);
mesh.y = [temp1(1:end-1), mesh.y, temp2(2:end)]; % add coarser discretization
mesh.z = linspace(0,substrate_thickness,5); % discretize the substrate with 5 cells
temp1 = linspace(substrate_thickness,2*substrate_thickness,5);
mesh.z = [mesh.z temp1(2:end)]; % add same space above the strip
temp1 = linspace(2*substrate_thickness,5*substrate_thickness,10);
mesh.z = [mesh.z temp1(2:end)]; % coarser discretization
CSX = DefineRectGrid( CSX, drawingunits, mesh );

%% Material definitions
CSX = AddMetal( CSX, 'PEC' );
CSX = AddMaterial( CSX, 'RO4350B' );

%% substrate
CSX = SetMaterialProperty( CSX, 'RO4350B', 'Epsilon', 3.66 );
start = [mesh.x(1),   mesh.y(1),   0];
stop  = [mesh.x(end), mesh.y(end), substrate_thickness];
CSX = AddBox( CSX, 'RO4350B', 0, start, stop );

%% MSL port
CSX = AddExcitation( CSX, 'excite', 0, [0 0 1]);
portstart = [          0, -MSL_width/2, substrate_thickness];
portstop =  [ MSL_length,  MSL_width/2, 0];
[CSX,portstruct] = AddMSLPort( CSX, 1, 'PEC', portstart, portstop, [1 0 0], [0 0 1], 'excite' );

%% MSL
start = [-MSL_length, -MSL_width/2, substrate_thickness];
stop  = [          0,  MSL_width/2, substrate_thickness];
CSX = AddBox( CSX, 'PEC', 0, start, stop );

%% define dump boxes... %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CSX = AddDump(CSX,'Et_','DumpType',0,'DumpMode',0);
start = [mesh.x(1)  , mesh.y(1),   substrate_thickness/2];
stop  = [mesh.x(end), mesh.y(end), substrate_thickness/2];
CSX = AddBox(CSX,'Et_',0 , start,stop);

CSX = AddDump(CSX,'Ht_','DumpType',1,'DumpMode',0);
CSX = AddBox(CSX,'Ht_',0,start,stop);
 

%% define openEMS options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
openEMS_opts = '';
openEMS_opts = [openEMS_opts ' --disable-dumps'];
% openEMS_opts = [openEMS_opts ' --debug-material'];
% openEMS_opts = [openEMS_opts ' --debug-operator'];
% openEMS_opts = [openEMS_opts ' --debug-boxes'];
% openEMS_opts = [openEMS_opts ' --engine=sse-compressed'];
% openEMS_opts = [openEMS_opts ' --engine=multithreaded'];
openEMS_opts = [openEMS_opts ' --engine=fastest'];

Sim_Path = 'tmp';
Sim_CSX = 'MSL2.xml';

if ~postprocessing_only
    rmdir(Sim_Path,'s');
end
mkdir(Sim_Path);

%% Write openEMS compatible xml-file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WriteOpenEMS( [Sim_Path '/' Sim_CSX], FDTD, CSX );

%% cd to working dir and run openEMS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
savePath = pwd;
cd(Sim_Path); %cd to working dir
args = [Sim_CSX ' ' openEMS_opts];
if ~postprocessing_only
    invoke_openEMS(args);
end
cd(savePath);



%% postproc & do the plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U = ReadUI({'port_ut1A','port_ut1B','et'},'tmp/');
I = ReadUI({'port_it1A','port_it1B'},'tmp/');
delta_t_2 = I.TD{1}.t(1) - U.TD{1}.t(1); % half time-step (s) 

% create finer frequency resolution
f = linspace( 0, f_max, 1601 );
for n=1:numel(U.FD)
	U.FD{n}.f = f;
	U.FD{n}.val = DFT_time2freq( U.TD{n}.t, U.TD{n}.val, f );
end
for n=1:numel(I.FD)
	I.FD{n}.f = f;
	I.FD{n}.val = DFT_time2freq( I.TD{n}.t, I.TD{n}.val, f );
	I.FD{n}.val = I.FD{n}.val .* exp(-1i*2*pi*I.FD{n}.f*delta_t_2); % compensate half time-step advance of H-field
end

% interpolate et to the time spacing of the voltage probes
et = interp1( U.TD{3}.t, U.TD{3}.val, U.TD{1}.t );

f = U.FD{1}.f;

% Z = (U.FD{1}.val+U.FD{2}.val)/2 ./ I.FD{1}.val;
% plot( f*1e-9, [real(Z);imag(Z)],'Linewidth',2);
% xlabel('frequency (GHz)');
% ylabel('impedance (Ohm)');
% grid on;
% legend( {'real','imaginary'}, 'location', 'northwest' )
% title( 'line impedance (will fail in case of reflections!)' );

% figure
% plotyy(U.TD{1}.t/1e-6,[U.TD{1}.val;U.TD{2}.val],U.TD{1}.t/1e-6,et);
% xlabel('time (us)');
% ylabel('amplitude (V)');
% grid on;
% title( 'Time domain voltage probes and excitation signal' );
% 
% figure
% plot(I.TD{1}.t/1e-6,[I.TD{1}.val;I.TD{2}.val]);
% xlabel('time (us)');
% ylabel('amplitude (A)');
% grid on;
% title( 'Time domain current probes' );


%% port analysis
[S11,beta,ZL] = calcMSLPort( portstruct, Sim_Path, f );

figure
plot( sin(0:0.01:2*pi), cos(0:0.01:2*pi), 'Color', [.7 .7 .7] );
hold on
plot( 0.5+0.5*sin(0:0.01:2*pi), 0.5*cos(0:0.01:2*pi), 'Color', [.7 .7 .7] );
plot( [-1 1], [0 0], 'Color', [.7 .7 .7] );
plot( S11, 'k' );
plot( real(S11(1)), imag(S11(1)), '*r' );
axis equal
title( 'Reflection coefficient S11 at the measurement plane' );

figure
plot( f/1e9, [real(S11);imag(S11)], 'Linewidth',2 );
legend( {'Re(S11)', 'Im(S11)'} );
ylabel( 'amplitude' );
xlabel( 'frequency (GHz)' );
title( 'Reflection coefficient S11 at the measurement plane' );

figure
plotyy( f/1e9, 20*log10(abs(S11)), f/1e9, angle(S11)/pi*180 );
legend( {'abs(S11)', 'angle(S11)'} );
xlabel( 'frequency (GHz)' );
title( 'Reflection coefficient S11 at the measurement plane' );

figure
plot( f/1e9, [real(beta);imag(beta)], 'Linewidth',2 );
legend( 'Re(beta)', 'Im(beta)' );
ylabel( 'propagation constant beta (1/m)' );
xlabel( 'frequency (GHz)' );
title( 'Propagation constant of the MSL' );

figure
plot( f/1e9, [real(ZL);imag(ZL)], 'Linewidth',2);
xlabel('frequency (GHz)');
ylabel('impedance (Ohm)');
grid on;
legend( {'real','imaginary'}, 'location', 'northeast' )
title( 'Characteristic line impedance ZL' );
ylim( [-2*mean(real(ZL)) 2*mean(real(ZL))] );

% reference plane shift (to the end of the port)
ref_shift = abs(portstop(1) - portstart(1));
[S11,beta,ZL] = calcMSLPort( portstruct, Sim_Path, f, ref_shift );

figure
plotyy( f/1e9, 20*log10(abs(S11)), f/1e9, angle(S11)/pi*180 );
legend( {'abs(S11)', 'angle(S11)'} );
xlabel( 'frequency (GHz)' );
title( 'Reflection coefficient S11 at the reference plane (at the electric wall)' );

figure
plot( sin(0:0.01:2*pi), cos(0:0.01:2*pi), 'Color', [.7 .7 .7] );
hold on
plot( 0.5+0.5*sin(0:0.01:2*pi), 0.5*cos(0:0.01:2*pi), 'Color', [.7 .7 .7] );
plot( [-1 1], [0 0], 'Color', [.7 .7 .7] );
plot( S11, 'k' );
plot( real(S11(1)), imag(S11(1)), '*r' );
axis equal
title( 'Reflection coefficient S11 at the reference plane (at the electric wall)' );
