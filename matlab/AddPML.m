function mesh = AddPML( mesh, numcells, CoordSystem )
% mesh = AddPML( mesh, numcells, <CoordSystem>  )
%
% Adds equidistant cells to the specified directions of the simulation
% area. This is used to put a PML (perfectly matched layer) absorber there.
% Remember: this function only adds space for the PML, the boundary
% conditions need to be set correctly to really add PML material.
%
% The mesh is sorted and duplicate lines are removed.
%
% input:
%   mesh:     mesh structure
%   numcells: 1x6 vector (xmin,xmax,ymin,ymax,zmin,zmax) with number of
%             cells to add to this direction
%   CoordSystem (optional): set to 1 in case of cylindrical mesh using
%               mesh.r, mesh.a and mesh.z
%
% output:
%   mesh:     new mesh with the added lines
%
% example:
%   % some fixed mesh lines
%   mesh.x = [-100 0 100];
%   mesh.y = [-100 0 100];
%   mesh.z = [0 500];
%   mesh = DetectEdges(CSX, mesh); %detect edges
%   mesh = SmoothMesh(mesh, c0/fmax/20/unit); % smooth the mesh
%
%   mesh = AddPML(mesh, 8); % add 8 lines to all directions
%   % or
%   mesh = AddPML(mesh, [0 0 0 0 8 8]); % add 8 lines in both z-directions
%
% See also DefineRectGrid, SmoothMesh, DetectEdges
%
% openEMS matlab interface
% -----------------------
% Sebastian Held <sebastian.held@uni-due.de>

% check
error( nargchk(2,3,nargin) );

if (numel(numcells)==1)
    numcells = ones(6,1)*numcells;
end

numcells = reshape( numcells, 1, [] );
if numel(numcells) ~= 6
    error( 'argument numcells needs to have exactly 6 elements' );
end

if (nargin<3)
    CoordSystem = 0;
end

dir_names = 'xyz';
if (CoordSystem==1)
    dir_names = 'raz';
end

mesh.(dir_names(1)) = unique(sort(mesh.(dir_names(1))));
mesh.(dir_names(2)) = unique(sort(mesh.(dir_names(2))));
mesh.(dir_names(3)) = unique(sort(mesh.(dir_names(3))));

% xmin
delta  = mesh.(dir_names(1))(2) - mesh.(dir_names(1))(1);
start  = mesh.(dir_names(1))(1) - numcells(1)*delta;
mesh.(dir_names(1)) = [start:delta:(mesh.(dir_names(1))(1)-delta), mesh.(dir_names(1))];

% xmax
delta  = mesh.(dir_names(1))(end) - mesh.(dir_names(1))(end-1);
stop   = mesh.(dir_names(1))(end) + numcells(2)*delta;
mesh.(dir_names(1)) = [mesh.(dir_names(1)), (mesh.(dir_names(1))(end)+delta):delta:stop];

% ymin
delta  = mesh.(dir_names(2))(2) - mesh.(dir_names(2))(1);
start  = mesh.(dir_names(2))(1) - numcells(3)*delta;
mesh.(dir_names(2)) = [start:delta:(mesh.(dir_names(2))(1)-delta), mesh.(dir_names(2))];

% ymax
delta  = mesh.(dir_names(2))(end) - mesh.(dir_names(2))(end-1);
stop   = mesh.(dir_names(2))(end) + numcells(4)*delta;
mesh.(dir_names(2)) = [mesh.(dir_names(2)), (mesh.(dir_names(2))(end)+delta):delta:stop];

% zmin
delta  = mesh.(dir_names(3))(2) - mesh.(dir_names(3))(1);
start  = mesh.(dir_names(3))(1) - numcells(5)*delta;
mesh.(dir_names(3)) = [start:delta:(mesh.(dir_names(3))(1)-delta), mesh.(dir_names(3))];

% zmax
delta  = mesh.(dir_names(3))(end) - mesh.(dir_names(3))(end-1);
stop   = mesh.(dir_names(3))(end) + numcells(6)*delta;
mesh.(dir_names(3)) = [mesh.(dir_names(3)), (mesh.(dir_names(3))(end)+delta):delta:stop];
