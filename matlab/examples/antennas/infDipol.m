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
% enlarge the box to be sure that one mesh line is covered by it
start = start - [0.1 0.1 0.1] * dipole_length/2;
stop  = stop  + [0.1 0.1 0.1] * dipole_length/2;
CSX = AddBox( CSX, 'infDipole', 1, start, stop );

% NFFF contour
start = [mesh.x(1)   mesh.y(1)   mesh.z(1) ];
stop  = [mesh.x(end) mesh.y(end) mesh.z(end) ];
[CSX nf2ff] = CreateNF2FFBox(CSX, 'nf2ff', start, stop);

% add space for PML
mesh = AddPML( mesh, [8 8 8 8 8 8] );
% define the mesh
CSX = DefineRectGrid( CSX, drawingunit, mesh );

if ~postprocessing_only
    % setup FDTD parameters & excitation function
    max_timesteps = 2000;
    min_decrement = 1e-6;
    FDTD = InitFDTD( 'NrTS', max_timesteps, 'EndCriteria', min_decrement, 'OverSampling',10 );
    FDTD = SetGaussExcite( FDTD, f_max/2, f_max/2 );
    BC = {'PML_8' 'PML_8' 'PML_8' 'PML_8' 'PML_8' 'PML_8'};
    FDTD = SetBoundaryCond( FDTD, BC );

    % Write openEMS compatible xml-file
    [~,~,~] = rmdir(Sim_Path,'s');
    [~,~,~] = mkdir(Sim_Path);
    WriteOpenEMS([Sim_Path '/' Sim_CSX],FDTD,CSX);

    % take a view at the "structure"
    CSXGeomPlot( [Sim_Path '/' Sim_CSX] );

    % define openEMS options and start simulation
    openEMS_opts = '';
    RunOpenEMS( Sim_Path, Sim_CSX, openEMS_opts );
end

%% post processing
disp( ' ' );
disp( ' ********************************************************** ' );
disp( ' ' );

% calculate the far field at phi=0 degrees and at phi=90 degrees
thetaRange = 0:0.5:359;
disp( 'calculating far field at phi=[0 90] deg..' );
nf2ff = CalcNF2FF( nf2ff, Sim_Path, f_max, thetaRange/180*pi, [0 pi/2], 'Mode', 1 );
Prad = nf2ff.Prad;
Dmax = nf2ff.Dmax;

theta_HPBW = interp1(nf2ff.E_norm{1}(find(thetaRange<90),1)/max(nf2ff.E_norm{1}(find(thetaRange<90),1)),thetaRange(find(thetaRange<90)),1/sqrt(2))*2;

% display power and directivity
disp( ['radiated power: Prad = ' num2str(Prad)] );
disp( ['directivity: Dmax = ' num2str(Dmax)] );
disp( ['theta_HPBW = ' num2str(theta_HPBW) ' °']);

% display polar plot for the e-field magnitude for phi = 0 & 90 deg
figure
polarFF(nf2ff,'xaxis','theta','param',[1 2]);

%% calculate the far field at theta=90 degrees
phiRange = 0:2:359;
disp( 'calculating far field at theta=90 deg..' );
nf2ff = CalcNF2FF( nf2ff, Sim_Path, f_max, 90/180*pi, phiRange/180*pi, 'Mode', 1 );

% display polar plot
figure
polarFF(nf2ff,'xaxis','phi','param',1);

%% calculate 3D pattern
phiRange = 0:5:360;
thetaRange = 0:5:180;
disp( 'calculating 3D far field...' );
nf2ff = CalcNF2FF( nf2ff, Sim_Path, f_max, thetaRange/180*pi, phiRange/180*pi, 'Mode', 1 );
figure
plotFF3D(nf2ff)

%%
E_far_normalized = nf2ff.E_norm{1} / max(nf2ff.E_norm{1}(:));
DumpFF2VTK([Sim_Path '/FF_pattern.vtk'],E_far_normalized, thetaRange,  phiRange);
disp(['view the farfield pattern "' Sim_Path '/FF_pattern.vtk" using paraview' ]);
