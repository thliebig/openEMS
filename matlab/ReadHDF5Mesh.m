function hdf_mesh = ReadHDF5Mesh(file)
% function hdf_mesh = ReadHDF5Mesh(file)
%
% returns:
% hdf_mesh.type
% hdf_mesh.names
% hdf_mesh.lines
%
% openEMS matlab interface
% -----------------------
% author: Thorsten Liebig

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
    hdf_mesh.lines{n} = hdf5read(file,names{n});
end

if (strcmp(names{1},'/mesh/rho'))
    hdf_mesh.type=1;
else
    hdf_mesh.type=0;
end