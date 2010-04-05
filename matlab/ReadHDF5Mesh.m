function hdf_mesh = ReadHDF5Mesh(file)

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