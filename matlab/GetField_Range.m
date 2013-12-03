function [field_i mesh_i] = GetField_Range(field, mesh, range)
% [field_i mesh_i] = GetField_Range(field, mesh, range)
%
%   Get a field dump subset within a given mesh range 
% 
%   example:
%       % specify a mesh range
%       range{1} = [0 150] * 1e-3;  % x in range 0..150mm
%       range{2} = [0];             % only one line close to y==0
%       range{3} = [];              % no range restriction
% 
%       % read hdf data
%       [field mesh] = ReadHDF5Dump('Et.h5');
%       % extract a ranged subset
%       [field_i mesh_i] = GetField_Range(field, mesh, range);
%     
%       %or both steps in one with the same result:
%       [field_i mesh_i] = ReadHDF5Dump('Et.h5','Range', range);
% 
% openEMS matlab interface
% -----------------------
% author: Thorsten Liebig
%
% See also ReadHDF5Dump ReadHDF5FieldData ReadHDF5Mesh

mesh_i = mesh;
for n=1:3
    if (numel(range{n})==0)
        ind_range{n} = [];
        
        ind_range{n} = 1:numel( mesh.lines{n});

    elseif (numel(range{n})==1)
        ind_range{n} = find( mesh.lines{n}>=range{n}(1) , 1);
        
        if (isempty(ind_range{n}))
            ind_range{n} = find( mesh.lines{n}>=range{n}(1) , 1, 'first');
        end
        if (isempty(ind_range{n}))
            ind_range{n} = find( mesh.lines{n}<=range{n}(2) , 1, 'last');
        end

    else
        ind_range{n} = find( mesh.lines{n}>=range{n}(1) & mesh.lines{n}<=range{n}(2));
    end
    
    mesh_i.lines{n} = mesh.lines{n}(ind_range{n});
end

% store original indices
if (isfield(mesh_i,'original_indices'))
    for n=1:3
        mesh_i.original_indices{n} = mesh_i.original_indices{n}(ind_range{n});
    end
else
    mesh_i.original_indices = ind_range;
end

field_i = field;

if (isfield(field,'FD'))
    for n=1:numel(field.FD.values)
        field_i.FD.values{n} = field.FD.values{n}(ind_range{1},ind_range{2},ind_range{3},:);
    end
end

if (isfield(field,'TD'))
	for n=1:numel(field.TD.values)
        field_i.TD.values{n} = field.TD.values{n}(ind_range{1},ind_range{2},ind_range{3},:);
    end
end
