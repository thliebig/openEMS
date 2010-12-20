function infDipol
%
% infinitesimal dipole example
%

close all
clear
clc
drawnow


postprocessing_only = 0;




global g
setup

dipole_orientation = 3; % 1,2,3: x,y,z
CSX = createStructure(dipole_orientation);
if ~postprocessing_only
    writeCSX( CSX );
    CSXGeomPlot( [g.Sim_Path '/' g.Sim_CSX] )
    run;
end
postprocess;






function setup
global g
physical_constants

% setup the simulation
g.drawingunit = 1e-6; % specify everything in um
g.Sim_Path = 'tmp';
g.Sim_CSX = 'tmp.xml';

g.f_max  = 1e9;
g.lambda = c0/g.f_max;

% setup geometry values
g.dipole_length = g.lambda/50 /g.drawingunit;



function CSX = createStructure(dipole_orientation)
global g
physical_constants

CSX = InitCSX();

% create an equidistant mesh
mesh.x = -g.dipole_length*10:g.dipole_length/2:g.dipole_length*10;
mesh.y = -g.dipole_length*10:g.dipole_length/2:g.dipole_length*10;
mesh.z = -g.dipole_length*10:g.dipole_length/2:g.dipole_length*10;
mesh = AddPML( mesh, [8 8 8 8 8 8] ); % add space for PML
CSX = DefineRectGrid( CSX, g.drawingunit, mesh );

% excitation
ex_vector = [0 0 0];
ex_vector(dipole_orientation) = 1;
start = ex_vector * -g.dipole_length/2;
stop  = ex_vector *  g.dipole_length/2;
CSX = AddExcitation( CSX, 'infDipole', 1, ex_vector );
CSX = AddBox( CSX, 'infDipole', 1, start, stop );

% NFFF contour
s1 = [-4.5, -4.5, -4.5] * g.dipole_length/2;
s2 = [ 4.5,  4.5,  4.5] * g.dipole_length/2;
CSX = AddBox( AddDump(CSX,'Et_xn','DumpType',0,'DumpMode',2,'FileType',1), 'Et_xn', 0, s1, [s1(1) s2(2) s2(3)] );
CSX = AddBox( AddDump(CSX,'Et_xp','DumpType',0,'DumpMode',2,'FileType',1), 'Et_xp', 0, [s2(1) s1(2) s1(3)], s2 );
CSX = AddBox( AddDump(CSX,'Et_yn','DumpType',0,'DumpMode',2,'FileType',1), 'Et_yn', 0, s1, [s2(1) s1(2) s2(3)] );
CSX = AddBox( AddDump(CSX,'Et_yp','DumpType',0,'DumpMode',2,'FileType',1), 'Et_yp', 0, [s1(1) s2(2) s1(3)], s2 );
CSX = AddBox( AddDump(CSX,'Et_zn','DumpType',0,'DumpMode',2,'FileType',1), 'Et_zn', 0, s1, [s2(1) s2(2) s1(3)] );
CSX = AddBox( AddDump(CSX,'Et_zp','DumpType',0,'DumpMode',2,'FileType',1), 'Et_zp', 0, [s1(1) s1(2) s2(3)], s2 );
CSX = AddBox( AddDump(CSX,'Ht_xn','DumpType',1,'DumpMode',2,'FileType',1), 'Ht_xn', 0, s1, [s1(1) s2(2) s2(3)] );
CSX = AddBox( AddDump(CSX,'Ht_xp','DumpType',1,'DumpMode',2,'FileType',1), 'Ht_xp', 0, [s2(1) s1(2) s1(3)], s2 );
CSX = AddBox( AddDump(CSX,'Ht_yn','DumpType',1,'DumpMode',2,'FileType',1), 'Ht_yn', 0, s1, [s2(1) s1(2) s2(3)] );
CSX = AddBox( AddDump(CSX,'Ht_yp','DumpType',1,'DumpMode',2,'FileType',1), 'Ht_yp', 0, [s1(1) s2(2) s1(3)], s2 );
CSX = AddBox( AddDump(CSX,'Ht_zn','DumpType',1,'DumpMode',2,'FileType',1), 'Ht_zn', 0, s1, [s2(1) s2(2) s1(3)] );
CSX = AddBox( AddDump(CSX,'Ht_zp','DumpType',1,'DumpMode',2,'FileType',1), 'Ht_zp', 0, [s1(1) s1(2) s2(3)], s2 );






function writeCSX(CSX)
global g
% setup FDTD parameters & excitation function
max_timesteps = 2000;
min_decrement = 1e-6;
FDTD = InitFDTD( max_timesteps, min_decrement, 'OverSampling',10 );
FDTD = SetGaussExcite( FDTD, g.f_max/2, g.f_max/2 );
BC = {'PML_8' 'PML_8' 'PML_8' 'PML_8' 'PML_8' 'PML_8'};
FDTD = SetBoundaryCond( FDTD, BC );

% Write openEMS compatible xml-file
[~,~,~] = rmdir(g.Sim_Path,'s');
[~,~,~] = mkdir(g.Sim_Path);
WriteOpenEMS([g.Sim_Path '/' g.Sim_CSX],FDTD,CSX);




function run
global g
% define openEMS options and start simulation
openEMS_opts = '';
openEMS_opts = [openEMS_opts ' --engine=fastest'];
RunOpenEMS( g.Sim_Path, g.Sim_CSX, openEMS_opts );




function postprocess
global g

disp( ' ' );
disp( ' ********************************************************** ' );
disp( ' ' );

% NFFF contour
filenames_E = {'Et_xn.h5','Et_xp.h5','Et_yn.h5','Et_yp.h5','Et_zn.h5','Et_zp.h5'};
filenames_H = {'Ht_xn.h5','Ht_xp.h5','Ht_yn.h5','Ht_yp.h5','Ht_zn.h5','Ht_zp.h5'};

% calculate the far field at phi=0 degrees and at phi=90 degrees
thetaRange = 0:2:359;
r = 1; % evaluate fields at radius r
disp( 'calculating far field at phi=[0 90] deg...' );
[E_far_theta,E_far_phi,Prad,Dmax] = AnalyzeNFFF2( g.Sim_Path, filenames_E, filenames_H, g.f_max, thetaRange, [0 90], r );


% display power and directivity
disp( ['radiated power: Prad = ' num2str(Prad)] );
disp( ['directivity: Dmax = ' num2str(Dmax)] );


% calculate the e-field magnitude for phi = 0 deg
E_phi0_far = zeros(1,numel(thetaRange));
for n=1:numel(thetaRange)
    E_phi0_far(n) = norm( [E_far_theta(n,1) E_far_phi(n,1)] );
end

% display polar plot
figure
polar( thetaRange/180*pi, E_phi0_far );
ylabel( 'theta / deg' );
title( ['electrical far field (V/m) @r=' num2str(r) ' m  phi=0 deg'] );
legend( 'e-field magnitude', 'Location', 'BestOutside' );


% calculate the e-field magnitude for phi = 90 deg
E_phi90_far = zeros(1,numel(thetaRange));
for n=1:numel(thetaRange)
    E_phi90_far(n) = norm([E_far_theta(n,2) E_far_phi(n,2)]);
end

% display polar plot
figure
polar( thetaRange/180*pi, E_phi90_far );
ylabel( 'theta / deg' );
title( ['electrical far field (V/m) @r=' num2str(r) ' m  phi=90 deg'] );
legend( 'e-field magnitude', 'Location', 'BestOutside' );




% calculate the far field at theta=90 degrees
phiRange = 0:2:359;
r = 1; % evaluate fields at radius r
disp( 'calculating far field at theta=90 deg...' );
[E_far_theta,E_far_phi] = AnalyzeNFFF2( g.Sim_Path, filenames_E, filenames_H, g.f_max, 90, phiRange, r );

E_theta90_far = zeros(1,numel(phiRange));
for n=1:numel(phiRange)
    E_theta90_far(n) = norm([E_far_theta(1,n) E_far_phi(1,n)]);
end

% display polar plot
figure
polar( phiRange/180*pi, E_theta90_far );
ylabel( 'phi / deg' );
title( ['electrical far field (V/m) @r=' num2str(r) ' m  theta=90 deg'] );
legend( 'e-field magnitude', 'Location', 'BestOutside' );




% calculate 3D pattern
phiRange = 0:15:360;
thetaRange = 0:10:180;
r = 1; % evaluate fields at radius r
disp( 'calculating 3D far field...' );
[E_far_theta,E_far_phi] = AnalyzeNFFF2( g.Sim_Path, filenames_E, filenames_H, g.f_max, thetaRange, phiRange, r );
E_far = sqrt( abs(E_far_theta).^2 + abs(E_far_phi).^2 );
E_far_normalized = E_far / max(E_far(:));
[theta,phi] = ndgrid(thetaRange/180*pi,phiRange/180*pi);
x = E_far_normalized .* sin(theta) .* cos(phi);
y = E_far_normalized .* sin(theta) .* sin(phi);
z = E_far_normalized .* cos(theta);
figure
surf( x,y,z, E_far_normalized );
axis equal
xlabel( 'x' );
ylabel( 'y' );
zlabel( 'z' );
