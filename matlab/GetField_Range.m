function [field_i mesh_i] = GetField_Range(field, mesh, range)
% [field_i mesh_i] = GetField_Range(field, mesh, range)
%
%   Get a field dump subset within a given mesh range 
% 
%   example:
%       % specify a mesh range
%       range{1} = [0 150] * 1e-3;
%       range{2} = [0 0];
%       range{3} = [-850 -400] * 1e-3;
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
    ind_range{n} = find( mesh.lines{n}>=range{n}(1) & mesh.lines{n}<=range{n}(2));
    
    if (isempty(ind_range{n}))
        ind_range{n} = find( mesh.lines{n}>=range{n}(1) , 1, 'first');
    end
    if (isempty(ind_range{n}))
        ind_range{n} = find( mesh.lines{n}<=range{n}(2) , 1, 'last');
    end
    
    mesh_i.lines{n} = mesh.lines{n}(ind_range{n});
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
