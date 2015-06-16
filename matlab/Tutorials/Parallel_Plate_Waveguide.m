%
% Tutorials / Parallel_Plate_Waveguide
%
% Describtion at:
% http://openems.de/index.php/Tutorial:_Parallel_Plate_Waveguide
%
% Tested with
%  - Matlab 2011a / Octave 4.0
%  - openEMS v0.0.33
%
% (C) 2011,2012 Sebastian Held <sebastian.held@gmx.de>
% (C) 2011-2015 Thorsten Liebig <thorsten.liebig@gmx.de>

close all
clear
clc

% init and define FDTD parameter
FDTD = InitFDTD(100,0,'OverSampling',50);
FDTD = SetSinusExcite(FDTD,10e6);
BC = {'PMC' 'PMC' 'PEC' 'PEC' 'MUR' 'MUR'};
FDTD = SetBoundaryCond(FDTD,BC);

% init and define FDTD mesh
CSX = InitCSX();
mesh.x = -10:10;
mesh.y = -10:10;
mesh.z = -10:30;
CSX = DefineRectGrid(CSX, 1, mesh);

% define the excitation
CSX = AddExcitation(CSX,'excitation',0,[0 1 0]);
CSX = AddBox(CSX,'excitation',0,[-10 -10 0],[10 10 0]);

% define a time domain e-field dump box
CSX = AddDump(CSX,'Et','DumpMode',0);
CSX = AddBox(CSX,'Et',0,[-10 0 -10],[10 0 30]);

% remove old simulation results (if exist)
rmdir('tmp','s');mkdir('tmp');

% write openEMS xml data file
WriteOpenEMS('tmp/tmp.xml',FDTD,CSX);

% view defined structure
CSXGeomPlot( 'tmp/tmp.xml' );

% run openEMS simulation
RunOpenEMS('tmp','tmp.xml','');

disp('use Paraview to visualize the FDTD result...');
