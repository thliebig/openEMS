%
% EXAMPLE / other / lumped elements
%
% This example demonstrates how to:
%  - use lumped elements
% 
%
% Tested with
%  - Matlab 2009b
%  - openEMS v0.0.21-3
%
% (C) 2010 Thorsten Liebig <thorsten.liebig@uni-due.de>

close all
clear
clc

%% setup the simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_max = 100e6;
f_excite = 300e6;
SimBox = 100;
mesh_size = 2;

Lumped.R = 1000;
Lumped.C = 10e-12;

% the parasitice inductance of the feeding has to be deduced with a R=0
% simulation
parasitic_L = 63e-9;

%% define openEMS options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
openEMS_opts = '';
% openEMS_opts = [openEMS_opts ' --debug-material'];
% openEMS_opts = [openEMS_opts ' --debug-boxes'];
% openEMS_opts = [openEMS_opts ' --debug-operator'];

Sim_Path = 'tmp';
Sim_CSX = 'lumped.xml';

[status, message, messageid] = rmdir(Sim_Path,'s');
[status,message,messageid] = mkdir(Sim_Path);

%% setup FDTD parameter & excitation function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FDTD = InitFDTD(30000,1e-6);
FDTD = SetGaussExcite(FDTD,f_excite/2,f_excite/2);
BC = [1 1 1 1 1 1];
FDTD = SetBoundaryCond(FDTD,BC);

%% setup CSXCAD geometry & mesh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CSX = InitCSX();
mesh.x = SmoothMeshLines([-SimBox/2,+SimBox/2],mesh_size);
mesh.y = SmoothMeshLines([-SimBox/2,+SimBox/2],mesh_size);
mesh.z = SmoothMeshLines([-SimBox/2,+SimBox/2],mesh_size);
CSX = DefineRectGrid(CSX, 1e-3,mesh);


%% create structure
% insert curve port
start = [ 10 -10 0]; 
stop  = [ 10  10 0];
CSX = AddCurvePort(CSX,0,1,100,start,stop,'excite');

% insert lumped element
CSX = AddLumpedElement( CSX, 'Capacitor', 1, 'C', Lumped.C, 'R', Lumped.R);
start = [ -14 -4 -4]; 
stop  = [ -6  4  4];
CSX = AddBox( CSX, 'Capacitor', 0, start, stop );

% insert feeding wire
CSX = AddMetal(CSX,'metal');
%first point
points(1,1) = -10;
points(2,1) = 4;
points(3,1) = 0;
%second point
points(1,2) = -10;
points(2,2) = 15;
points(3,2) = 0;
%3 point
points(1,end+1) = 10;
points(2,end) = 15;
points(3,end) = 0;
%4 point
points(1,end+1) = 10;
points(2,end) = 10;
points(3,end) = 0;
CSX = AddCurve(CSX,'metal', 10, points);

points(2,:) = -1*points(2,:);
CSX = AddCurve(CSX,'metal', 10, points);

%% Write openEMS compatoble xml-file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WriteOpenEMS([Sim_Path '/' Sim_CSX],FDTD,CSX);

% CSXGeomPlot([Sim_Path '/' Sim_CSX]);

%% run openEMS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RunOpenEMS(Sim_Path, Sim_CSX, openEMS_opts);

%% postproc & do the plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = linspace(1e6,f_max,1001);
w = 2*pi*f;
% read currents and voltages
U = ReadUI('port_ut1','tmp/',f);
I = ReadUI('port_it1','tmp/',f);

% calculate analytic impedance
if (Lumped.R>=0)
    Z_a = Lumped.R*(1-1i*w*Lumped.C*Lumped.R)./(1+(w*Lumped.C*Lumped.R).^2);
else
    Z_a = -1i./(w*Lumped.C);
end
   
% calculate numerical impedance
Z = U.FD{1}.val./I.FD{1}.val;

% remove parasitic feeding effects
Z = Z - 1i*w*parasitic_L;

L = imag(Z)./w;
C = -1./(w.*imag(Z));
C(find(C<0)) = nan;
L(find(L<0)) = nan;
R = real(Z);

subplot(2,1,1);
plot(f*1e-6,C*1e12,'Linewidth',2);
xlabel('frequency (MHz)');
ylabel('capacitance (pF)');
grid on;
subplot(2,1,2);
plot(f*1e-6,L*1e9,'Linewidth',2);
xlabel('frequency (MHz)');
ylabel('inductance (nH)');
grid on;

figure();
plot(f*1e-6,R,'Linewidth',2);
hold on
plot(f*1e-6,imag(Z),'r--','Linewidth',2);

plot(f*1e-6,real(Z_a),'g-.','Linewidth',1);
plot(f*1e-6,imag(Z_a),'m--','Linewidth',1);

xlabel('frequency (MHz)');
ylabel('resistance (Ohm)');
grid on;
legend( '\Re\{Z\}','\Im\{Z\}','\Re\{Z_{analytisch}\}','\Im\{Z_{analytisch}\}', 'location', 'northeast' )

figure();
errorR = (R-real(Z_a))./R*100;
errorX = (imag(Z)-imag(Z_a))./imag(Z)*100;
plot(f*1e-6,errorR,'Linewidth',2);
hold on
grid on;
plot(f*1e-6,errorX,'r--','Linewidth',2);
xlabel('frequency (MHz)');
ylabel('error (%)');
