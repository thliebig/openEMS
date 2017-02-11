%
% Tutorials / CylindricalWave_CC
%
% Describtion at:
% http://openems.de/index.php/Tutorial:_2D_Cylindrical_Wave
%
% Tested with
%  - Matlab 2011a/ Octave 4.0
%  - openEMS v0.0.33
%
% (C) 2011-2015 Thorsten Liebig <thorsten.liebig@gmx.de>

close all
clear
clc

%% setup the simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
physical_constants
mesh_res = 10;      %desired mesh resolution
radius = 2560;      %simulation domain radius
split = ['80,160,320,640,1280']; %radii to split the mesh into sub-grids
split_N = 5;        %number of nested sub-grids
heigth = mesh_res*4;

f0 = 1e9;

exite_offset = 1300;
excite_angle = 45;

%% setup FDTD parameter & excitation function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FDTD = InitFDTD('NrTS', 100000, 'EndCriteria', 1e-4, 'CoordSystem', 1, 'MultiGrid', split);
FDTD = SetGaussExcite(FDTD,f0,f0/2);
BC = [0 3 0 0 0 0];             % pml in positive r-direction
FDTD = SetBoundaryCond(FDTD,BC);

%% setup CSXCAD geometry & mesh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 50 mesh lines for the inner most mesh
% increase the total number of meshlines in alpha direcion for all sub-grids
N_alpha = 50 * 2^split_N + 1;

CSX = InitCSX('CoordSystem',1);
mesh.r = SmoothMeshLines([0 radius],mesh_res);
mesh.a = linspace(-pi,pi,N_alpha);
mesh.z = SmoothMeshLines([-heigth/2 0 heigth/2],mesh_res);
CSX = DefineRectGrid(CSX, 1e-3,mesh);

%% add the dipol %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
start = [exite_offset excite_angle/180*pi-0.001 -20];
stop =  [exite_offset excite_angle/180*pi+0.001  20];
if (exite_offset==0)
    start(2) = mesh.a(1);
    stop(2)  = mesh.a(1);
end
CSX = AddExcitation(CSX,'excite',1,[0 0 1]);
CSX = AddBox(CSX,'excite',0 ,start,stop);

%% define dump boxes... %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
start = [mesh.r(1)   mesh.a(1)   0];
stop =  [mesh.r(end-8) mesh.a(end) 0];

% time domain vtk dump
CSX = AddDump(CSX,'Et_ra','DumpType',0,'FileType',0,'SubSampling','4,10,1');
CSX = AddBox(CSX,'Et_ra',0 , start,stop);

% frequency domain hdf5 dump
CSX = AddDump(CSX,'Ef_ra','DumpType',10,'FileType',1,'SubSampling','2,2,2','Frequency',f0);
CSX = AddBox(CSX,'Ef_ra',0 , start,stop);

%% write/run the openEMS compatible xml-file
Sim_Path = 'tmp';
Sim_CSX = '2D_CC_Wave.xml';

[status, message, messageid] = rmdir( Sim_Path, 's' ); % clear previous directory
[status, message, messageid] = mkdir( Sim_Path ); % create empty simulation folder

WriteOpenEMS([Sim_Path '/' Sim_CSX],FDTD,CSX);
RunOpenEMS(Sim_Path, Sim_CSX);

%%
disp('use Paraview to visualize the vtk field dump...');

%%
[field mesh_h5] = ReadHDF5Dump([Sim_Path '/Ef_ra.h5']);

r = mesh_h5.lines{1};
a = mesh_h5.lines{2};
a(end+1) = a(1);            %closeup mesh for visualization
[R A] = ndgrid(r,a);
X = R.*cos(A);
Y = R.*sin(A);

Ez = squeeze(field.FD.values{1}(:,:,1,3));
Ez(:,end+1) = Ez(:,1);      %closeup mesh for visualization

E_max = max(max(abs(Ez)));  %get maximum E_z amplitude

while 1
    for ph = linspace(0,360,41) %animate phase from 0..360 degree
        surf(X,Y,real(Ez*exp(1j*ph*pi/180)),'EdgeColor','none')
        caxis([-E_max E_max]/10)
        zlim([-E_max E_max])
        pause(0.3)
    end
end
