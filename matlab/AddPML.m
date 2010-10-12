function mesh = AddPML( mesh, numcells )
%mesh = AddPML( mesh, numcells )
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
%     .x:     1xn vector (lines in x direction)
%     .y:     1xn vector (lines in y direction)
%     .z:     1xn vector (lines in z direction)
%   numcells: 1x6 vector (xmin,xmax,ymin,ymax,zmin,zmax) with number of
%             cells to add to this direction
%
% output:
%   mesh:     new mesh
%
% openEMS matlab interface
% -----------------------
% Sebastian Held <sebastian.held@uni-due.de>
% See also DefineRectGrid

% check
error( nargchk(2,2,nargin) );

numcells = reshape( numcells, 1, [] );
if numel(numcells) ~= 6
    error( 'argument numcells needs to have exactly 6 elements' );
end

mesh.x = unique(sort(mesh.x));
mesh.y = unique(sort(mesh.y));
mesh.z = unique(sort(mesh.z));

% xmin
delta  = mesh.x(2) - mesh.x(1);
start  = mesh.x(1) - numcells(1)*delta;
mesh.x = [start:delta:(mesh.x(1)-delta), mesh.x];

% xmax
delta  = mesh.x(end) - mesh.x(end-1);
stop   = mesh.x(end) + numcells(2)*delta;
mesh.x = [mesh.x, (mesh.x(end)+delta):delta:stop];

% ymin
delta  = mesh.y(2) - mesh.y(1);
start  = mesh.y(1) - numcells(3)*delta;
mesh.y = [start:delta:(mesh.y(1)-delta), mesh.y];

% ymax
delta  = mesh.y(end) - mesh.y(end-1);
stop   = mesh.y(end) + numcells(4)*delta;
mesh.y = [mesh.y, (mesh.y(end)+delta):delta:stop];

% zmin
delta  = mesh.z(2) - mesh.z(1);
start  = mesh.z(1) - numcells(5)*delta;
mesh.z = [start:delta:(mesh.z(1)-delta), mesh.z];

% zmax
delta  = mesh.z(end) - mesh.z(end-1);
stop   = mesh.z(end) + numcells(6)*delta;
mesh.z = [mesh.z, (mesh.z(end)+delta):delta:stop];
