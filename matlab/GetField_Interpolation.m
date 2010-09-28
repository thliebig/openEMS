function [field_i mesh_i] = GetField_Interpolation(field, mesh, numLines, varargin)
% [field_i mesh_i] = GetField_Interpolation(field, mesh, numLines, varargin)
%
%   example:
%   [field mesh] = ReadHDF5Dump('Et.h5');
%   %interpolate on a mesh with 21x21x101 lines
%   [field_i mesh_i] = GetField_Interpolation(field, mesh, [21 21 101]);
%   
%   %or both steps in one with the same result:
%   [field_i mesh_i] = ReadHDF5Dump('Et.h5','Interpolation', [21 21 101]);
%
% openEMS matlab interface
% -----------------------
% author: Thorsten Liebig
%
% See also ReadHDF5Dump ReadHDF5FieldData ReadHDF5Mesh

if (~isnumeric(numLines) || numel(numLines)~=3)
    error('openEMS:GetField_Interpolation: numLines for interpolation must be a vector...');
end
    
x = mesh.lines{1};
y = mesh.lines{2};
z = mesh.lines{3};

x_i = linspace(x(1),x(end),numLines(1));
y_i = linspace(y(1),y(end),numLines(2));
z_i = linspace(z(1),z(end),numLines(3));

field_i = field;
mesh_i = mesh;
mesh_i.lines{1} = x_i;
mesh_i.lines{2} = y_i;
mesh_i.lines{3} = z_i;

NULL = zeros(numel(x_i),numel(y_i),numel(z_i),3);
for n=1:numel(field.values)
    field_i.values{n} = NULL;
end
clear NULL;

% matlab cannot handle 3D data to be 2D data, workaround for these cases
if (numel(x)==1)
    [Y Z] = ndgrid(y,z);
    [Y_I Z_I] = ndgrid(y_i,z_i);
    for n=1:numel(field.values)
        field_i.values{n}(1,:,:,1) = interpn(Y,Z,squeeze(field.values{n}(1,:,:,1)),Y_I,Z_I);
        field_i.values{n}(1,:,:,2) = interpn(Y,Z,squeeze(field.values{n}(1,:,:,2)),Y_I,Z_I);
        field_i.values{n}(1,:,:,3) = interpn(Y,Z,squeeze(field.values{n}(1,:,:,3)),Y_I,Z_I);
    end
    return;
end

if (numel(y)==1)
    [X Z] = ndgrid(x,z);
    [X_I Z_I] = ndgrid(x_i,z_i);
    for n=1:numel(field.values)
        field_i.values{n}(:,1,:,1) = interpn(X,Z,squeeze(field.values{n}(:,1,:,1)),X_I,Z_I);
        field_i.values{n}(:,1,:,2) = interpn(X,Z,squeeze(field.values{n}(:,1,:,2)),X_I,Z_I);
        field_i.values{n}(:,1,:,3) = interpn(X,Z,squeeze(field.values{n}(:,1,:,3)),X_I,Z_I);
    end
    return;
end

if (numel(z)==1)
    [X Y] = ndgrid(x,y);
    [X_I Y_I] = ndgrid(x_i,y_i);
    for n=1:numel(field.values)
        field_i.values{n}(:,:,1,1) = interpn(X,Y,squeeze(field.values{n}(:,:,1,1)),X_I,Y_I);
        field_i.values{n}(:,:,1,2) = interpn(X,Y,squeeze(field.values{n}(:,:,1,2)),X_I,Y_I);
        field_i.values{n}(:,:,1,3) = interpn(X,Y,squeeze(field.values{n}(:,:,1,3)),X_I,Y_I);
    end
    return;
end


%real 3D case
[X Y Z] = ndgrid(x,y,z);
[X_I Y_I Z_I] = ndgrid(x_i,y_i,z_i);
for n=1:numel(field.values)
    field_i.values{n}(:,:,:,1) = interpn(X,Y,Z,field.values{n}(:,:,:,1),X_I,Y_I,Z_I);
    field_i.values{n}(:,:,:,2) = interpn(X,Y,Z,field.values{n}(:,:,:,2),X_I,Y_I,Z_I);
    field_i.values{n}(:,:,:,3) = interpn(X,Y,Z,field.values{n}(:,:,:,3),X_I,Y_I,Z_I);
end
