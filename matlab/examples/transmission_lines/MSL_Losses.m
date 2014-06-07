%
% examples / microstrip / MSL_Losses
%
% This example demonstrates how to model sheet conductor losses
%
% Tested with
%  - Matlab 2013a / Octave 3.8.1+
%  - openEMS v0.0.32
%
% (C) 2012-2014 Thorsten Liebig <thorsten.liebig@gmx.de>

close all
clear
clc

%% setup the simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
physical_constants;
unit = 1e-6; % specify everything in um
MSL.length = 10000;
MSL.port_dist = 5000;
MSL.width = 225;
MSL.conductivity = 41e6;
MSL.thickness = 35e-6;

substrate.thickness = 250;
substrate.epr = 9.8;

f_start = 0e9;
f_stop  = 25e9;

lambda = c0/f_stop;

%% setup FDTD parameters & excitation function %%%%%%%%%%%%%%%%%%%%%%%%%%%%
FDTD = InitFDTD('endCriteria',1e-4);
FDTD = SetGaussExcite(FDTD,0.5*(f_start+f_stop),0.5*(f_stop-f_start));
BC   = {'PML_8' 'PML_8' 'PML_8' 'PML_8' 'PEC' 'PML_8'};
FDTD = SetBoundaryCond( FDTD, BC );

%% setup CSXCAD geometry & mesh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CSX = InitCSX();
resolution = c0/(f_stop*sqrt(substrate.epr))/unit /20;
mesh.x = SmoothMeshLines( [-MSL.length*0.5-MSL.port_dist 0 MSL.length*0.5+MSL.port_dist], resolution, 1.3 ,0 );
mesh.y = SmoothMeshLines2( [0 MSL.width/2], resolution/6 , 1.3);
mesh.y = SmoothMeshLines( [-0.5*lambda/unit -mesh.y mesh.y 0.5*lambda/unit], resolution, 1.4);
mesh.z = SmoothMeshLines( [linspace(0,substrate.thickness,10) 0.5*lambda/unit], resolution );
CSX = DefineRectGrid( CSX, unit, mesh );

%% substrate
CSX = AddMaterial( CSX, 'RO4350B' );
CSX = SetMaterialProperty( CSX, 'RO4350B', 'Epsilon', substrate.epr );
start = [mesh.x(1),   mesh.y(1),   0];
stop  = [mesh.x(end), mesh.y(end), substrate.thickness];
CSX = AddBox( CSX, 'RO4350B', 0, start, stop );

%% MSL ports and lossy line
CSX = AddConductingSheet( CSX, 'gold', MSL.conductivity, MSL.thickness );
portstart = [ mesh.x(1),               -MSL.width/2, substrate.thickness];
portstop  = [ mesh.x(1)+MSL.port_dist,  MSL.width/2, 0];
[CSX, port{1}] = AddMSLPort( CSX, 999, 1, 'gold', portstart, portstop, 0, [0 0 -1], 'ExcitePort', true, 'FeedShift', 10*resolution, 'MeasPlaneShift',  MSL.port_dist);

portstart = [mesh.x(end),              -MSL.width/2, substrate.thickness];
portstop  = [mesh.x(end)-MSL.port_dist, MSL.width/2, 0];
[CSX, port{2}] = AddMSLPort( CSX, 999, 2, 'gold', portstart, portstop, 0, [0 0 -1], 'MeasPlaneShift',  MSL.port_dist );

start = [mesh.x(1)+MSL.port_dist,   -MSL.width/2, substrate.thickness];
stop  = [mesh.x(end)-MSL.port_dist,  MSL.width/2, substrate.thickness];
CSX = AddBox(CSX,'gold',500,start,stop);

%% write/show/run the openEMS compatible xml-file
Sim_Path = 'tmp';
Sim_CSX = 'msl.xml';

[status, message, messageid] = rmdir( Sim_Path, 's' ); % clear previous directory
[status, message, messageid] = mkdir( Sim_Path ); % create empty simulation folder

WriteOpenEMS( [Sim_Path '/' Sim_CSX], FDTD, CSX );
CSXGeomPlot( [Sim_Path '/' Sim_CSX] );
RunOpenEMS( Sim_Path, Sim_CSX ,'');

%% post-processing
close all
f = linspace( f_start, f_stop, 1601 );
port = calcPort(port, Sim_Path, f, 'RefImpedance', 50);

s11 = port{1}.uf.ref./ port{1}.uf.inc;
s21 = port{2}.uf.ref./ port{1}.uf.inc;

plot(f/1e9,-20*log10(abs(s21)),'r--','LineWidth',2);
grid on;
hold on;
ylabel('-|S_21| (dB)','Interpreter','None');
xlabel('frequency (GHz)');

%% plot 35um thickness loss model curve
% values extracted from http://wcalc.sourceforge.net/cgi-bin/microstrip.cgi
model.f    = [1   2   2.5 3   4   5   7.5 10   12.5  15   17.5 20   25   ]; % frequency in GHz
model.loss = [3.0 4.2 4.7 5.2 5.9 6.6 8.1 9.38 10.5  11.5 12.4 13.2 14.65]; % loss in db/m

plot(model.f, model.loss * MSL.length * unit ,'k-','LineWidth',1);
legend('FDTD simulated attenuation','t=35um, loss model by E. Hammerstad & F. Bekkadal','Location','NorthWest');


