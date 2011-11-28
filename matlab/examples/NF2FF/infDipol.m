%
% infinitesimal dipole example
%

close all
clear
clc

postprocessing_only = 0;

physical_constants

% setup the simulation
drawingunit = 1e-6; % specify everything in um
Sim_Path = 'tmp';
Sim_CSX = 'tmp.xml';

f_max  = 1e9;
lambda = c0/f_max;

% setup geometry values
dipole_length = lambda/50 /drawingunit;


dipole_orientation = 3; % 1,2,3: x,y,z


CSX = InitCSX();

% create an equidistant mesh
mesh.x = -dipole_length*10:dipole_length/2:dipole_length*10;
mesh.y = -dipole_length*10:dipole_length/2:dipole_length*10;
mesh.z = -dipole_length*10:dipole_length/2:dipole_length*10;

% excitation
ex_vector = [0 0 0];
ex_vector(dipole_orientation) = 1;
start = ex_vector * -dipole_length/2;
stop  = ex_vector *  dipole_length/2;
CSX = AddExcitation( CSX, 'infDipole', 1, ex_vector );
CSX = AddBox( CSX, 'infDipole', 1, start, stop );

% NFFF contour
start = [mesh.x(1)   mesh.y(1)   mesh.z(1) ]
stop  = [mesh.x(end) mesh.y(end) mesh.z(end) ]
[CSX nf2ff] = CreateNF2FFBox(CSX, 'nf2ff', start, stop);

% add space for PML
mesh = AddPML( mesh, [8 8 8 8 8 8] );
% define the mesh
CSX = DefineRectGrid( CSX, drawingunit, mesh );

if ~postprocessing_only
    % setup FDTD parameters & excitation function
    max_timesteps = 2000;
    min_decrement = 1e-6;
    FDTD = InitFDTD( max_timesteps, min_decrement, 'OverSampling',10 );
    FDTD = SetGaussExcite( FDTD, f_max/2, f_max/2 );
    BC = {'PML_8' 'PML_8' 'PML_8' 'PML_8' 'PML_8' 'PML_8'};
    FDTD = SetBoundaryCond( FDTD, BC );

    % Write openEMS compatible xml-file
    [~,~,~] = rmdir(Sim_Path,'s');
    [~,~,~] = mkdir(Sim_Path);
    WriteOpenEMS([Sim_Path '/' Sim_CSX],FDTD,CSX);

    % define openEMS options and start simulation
    openEMS_opts = '';
    RunOpenEMS( Sim_Path, Sim_CSX, openEMS_opts );
end

%% post processing
disp( ' ' );
disp( ' ********************************************************** ' );
disp( ' ' );

% calculate the far field at phi=0 degrees and at phi=90 degrees
thetaRange = 0:2:359;
r = 1; % evaluate fields at radius r
disp( 'calculating far field at phi=[0 90] deg..' );
[E_far_theta,E_far_phi,Prad,Dmax] = AnalyzeNF2FF( Sim_Path, nf2ff, f_max, thetaRange, [0 90], r );

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
disp( 'calculating far field at theta=90 deg..' );
[E_far_theta,E_far_phi] = AnalyzeNF2FF( Sim_Path, nf2ff, f_max, 90, phiRange, r );

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
[E_far_theta,E_far_phi] = AnalyzeNF2FF( Sim_Path, nf2ff, f_max, thetaRange, phiRange, r );
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

%
DumpFF2VTK([Sim_Path '/FF_pattern.vtk'],E_far_normalized, thetaRange,  phiRange);
disp(['view the farfield pattern "' Sim_Path '/FF_pattern.vtk" using paraview' ]);
