function [field_i mesh_i] = GetField_Interpolation(field, mesh, lines, varargin)
% [field_i mesh_i] = GetField_Interpolation(field, mesh, lines, varargin)
%
%   Get an interpolated field, e.g. read by ReadHDF5Dump 
% 
%   homogeneous interpolation given by a 3x1 vector: e.g. [21,1,101]
% 
%   arbitrary interpolation on a given mesh:
%               e.g.:   mesh_interp{1} = linspace(0,  1,101) * 1e-3;
%                       mesh_interp{2} = linspace(0,0.5, 51) * 1e-3;
%                       mesh_interp{3} = linspace(0,0.2, 21) * 1e-3;
% 
%   example:
%   [field mesh] = ReadHDF5Dump('Et.h5');
%   %interpolate on a mesh with 21x21x101 lines
%   [field_i mesh_i] = GetField_Interpolation(field, mesh, [21 21 101]);
%   or
%   [field_i mesh_i] = GetField_Interpolation(field, mesh, mesh_interp);
%   
%   %or both steps in one with the same result:
%   [field_i mesh_i] = ReadHDF5Dump('Et.h5','Interpolation', [21 21 101]);
%
% openEMS matlab interface
% -----------------------
% author: Thorsten Liebig
%
% See also ReadHDF5Dump ReadHDF5FieldData ReadHDF5Mesh

if ((~iscell(lines) && ~isnumeric(lines)) || numel(lines)~=3)
    error('openEMS:GetField_Interpolation: numLines for interpolation must be a vector...');
end
    
x = mesh.lines{1};
y = mesh.lines{2};
z = mesh.lines{3};

if (isnumeric(lines))
    if (lines(1)==0)
        x_i = x;
    else
        x_i = linspace(x(1),x(end),lines(1));
    end
    if (lines(2)==0)
        y_i = y;
    else
        y_i = linspace(y(1),y(end),lines(2));
    end
    if (lines(3)==0)
        z_i = z;
    else
        z_i = linspace(z(1),z(end),lines(3));
    end
else
    if isempty(lines{1})
        x_i = x;
    else
        x_i = lines{1};
    end
    if isempty(lines{2})
        y_i = y;
    else
        y_i = lines{2};
    end
    if isempty(lines{3})
        z_i = z;
    else
        z_i = lines{3};
    end
end

field_i = field;
mesh_i = mesh;
mesh_i.lines{1} = x_i;
mesh_i.lines{2} = y_i;
mesh_i.lines{3} = z_i;

% clear or create empty original indices list, since such do not make any
% sense with interpolated field values
mesh_i.original_indices = {};

if (isfield(field,'TD'))
    field_i.TD = interpolate_fields(field.TD,x,y,z, x_i, y_i, z_i);
    field_i.TD.time = field.TD.time;
    field_i.TD.names= field.TD.names;
end

if (isfield(field,'FD'))
    field_i.FD = interpolate_fields(field.FD,x,y,z, x_i, y_i, z_i);
    field_i.FD.frequency = field.FD.frequency;
    field_i.FD.DataType = field.FD.DataType;
end

return

function field_i = interpolate_fields(field, x,y,z, x_i, y_i, z_i)

% matlab cannot handle 3D data to be 2D data, workaround for these cases
if (numel(x)==1)
    [Y Z] = ndgrid(y,z);
    [Y_I Z_I] = ndgrid(y_i,z_i);
    for n=1:numel(field.values)
        field_i.values{n}(1,:,:,1) = interpn(Y,Z,squeeze(field.values{n}(1,:,:,1)),Y_I,Z_I);
        if (size(field.values{n},4)>1)
            field_i.values{n}(1,:,:,2) = interpn(Y,Z,squeeze(field.values{n}(1,:,:,2)),Y_I,Z_I);
            field_i.values{n}(1,:,:,3) = interpn(Y,Z,squeeze(field.values{n}(1,:,:,3)),Y_I,Z_I);
        end
    end
    return;
end

if (numel(y)==1)
    [X Z] = ndgrid(x,z);
    [X_I Z_I] = ndgrid(x_i,z_i);
    for n=1:numel(field.values)
        field_i.values{n}(:,1,:,1) = interpn(X,Z,squeeze(field.values{n}(:,1,:,1)),X_I,Z_I);
        if (size(field.values{n},4)>1)
            field_i.values{n}(:,1,:,2) = interpn(X,Z,squeeze(field.values{n}(:,1,:,2)),X_I,Z_I);
            field_i.values{n}(:,1,:,3) = interpn(X,Z,squeeze(field.values{n}(:,1,:,3)),X_I,Z_I);
        end
    end
    return;
end

if (numel(z)==1)
    [X Y] = ndgrid(x,y);
    [X_I Y_I] = ndgrid(x_i,y_i);
    for n=1:numel(field.values)
        field_i.values{n}(:,:,1,1) = interpn(X,Y,squeeze(field.values{n}(:,:,1,1)),X_I,Y_I);
        if (size(field.values{n},4)>1)
            field_i.values{n}(:,:,1,2) = interpn(X,Y,squeeze(field.values{n}(:,:,1,2)),X_I,Y_I);
            field_i.values{n}(:,:,1,3) = interpn(X,Y,squeeze(field.values{n}(:,:,1,3)),X_I,Y_I);
        end
    end
    return;
end


%real 3D case
[X Y Z] = ndgrid(x,y,z);
[X_I Y_I Z_I] = ndgrid(x_i,y_i,z_i);
for n=1:numel(field.values)
    field_i.values{n}(:,:,:,1) = interpn(X,Y,Z,field.values{n}(:,:,:,1),X_I,Y_I,Z_I);
    if (size(field.values{n},4)>1)
        field_i.values{n}(:,:,:,2) = interpn(X,Y,Z,field.values{n}(:,:,:,2),X_I,Y_I,Z_I);
        field_i.values{n}(:,:,:,3) = interpn(X,Y,Z,field.values{n}(:,:,:,3),X_I,Y_I,Z_I);
    end
end
