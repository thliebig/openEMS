function [field_i mesh_i] = GetField_SubSampling(field, mesh, subsampling, varargin)
% [field_i mesh_i] = GetField_SubSampling(field, mesh, subsampling, varargin)
%
%   Get a sub-sampled field, e.g. read by ReadHDF5Dump
%
%   sub-sampling e.g. skipping every second line in x/r direction: [2 1 1]
%
% openEMS matlab interface
% -----------------------
% author: Thorsten Liebig
%
% See also ReadHDF5Dump ReadHDF5FieldData ReadHDF5Mesh

if (~isnumeric(subsampling) || numel(subsampling)~=3)
    error('openEMS:GetField_Interpolation: numLines for interpolation must be a vector...');
end

x = mesh.lines{1};
y = mesh.lines{2};
z = mesh.lines{3};

ss_idx{1} = 1:subsampling(1):numel(x);
ss_idx{2} = 1:subsampling(2):numel(y);
ss_idx{3} = 1:subsampling(3):numel(z);

x_i = x(ss_idx{1});
y_i = y(ss_idx{2});
z_i = z(ss_idx{3});

field_i = field;
mesh_i = mesh;
mesh_i.lines{1} = x_i;
mesh_i.lines{2} = y_i;
mesh_i.lines{3} = z_i;

% store original indices
if (isfield(mesh_i,'original_indices'))
    for n=1:3
        mesh_i.original_indices{n} = mesh_i.original_indices{n}(ss_idx{n});
    end
else
    mesh_i.original_indices = ss_idx;
end

if (isfield(field,'TD'))
    field_i.TD = subsample_fields(field.TD,ss_idx);
    field_i.TD.time = field.TD.time;
    field_i.TD.names= field.TD.names;
end

if (isfield(field,'FD'))
    field_i.FD = subsample_fields(field.FD,ss_idx);
    field_i.FD.frequency = field.FD.frequency;
    field_i.FD.DataType = field.FD.DataType;
end

return

function field_i = subsample_fields(field, ss_idx)

for n=1:numel(field.values)
    field_i.values{n} = field.values{n}(ss_idx{1},ss_idx{2},ss_idx{3},:);
end
