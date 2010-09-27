function hdf_mesh = ReadHDF5Mesh(file)
% function hdf_mesh = ReadHDF5Mesh(file)
%
%   Get the raw mesh data stored in the hdf5 dump file created by openEMS
%
% returns:
% hdf_mesh.type     (0-> cartesian, 1-> cylindrical mesh type)
% hdf_mesh.names    (e.g. 'Mesh/y')
% hdf_mesh.lines    (e.g. [0,1,2,3,4])
%
% openEMS matlab interface
% -----------------------
% author: Thorsten Liebig
%
% See also ReadHDF5FieldData

info = hdf5info(file);

for n=1:numel(info.GroupHierarchy.Groups)
    if strcmp(info.GroupHierarchy.Groups(n).Name,'/Mesh')
        for m=1:numel(info.GroupHierarchy.Groups(n).Datasets)
            names{m} = info.GroupHierarchy.Groups(n).Datasets(m).Name;
        end
    end
end

hdf_mesh.names = names;
for n=1:numel(names)
    hdf_mesh.lines{n} = double(hdf5read(file,names{n}));
end

if (strcmp(names{1},'/Mesh/alpha'))
    % alpha and rho are in the wrong order, flip to have rho, alpha, z
    hdf_mesh.names(1:2) = fliplr(hdf_mesh.names(1:2));
    hdf_mesh.lines(1:2) = fliplr(hdf_mesh.lines(1:2));
    hdf_mesh.type=1;
    return
end

hdf_mesh.type=0;
