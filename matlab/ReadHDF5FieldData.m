function hdf_fielddata = ReadHDF5FieldData(file)
% function hdf_fielddata = ReadHDF5FieldData(file)
%
% returns:
% hdf_fielddata.names
% hdf_fielddata.values
%
% openEMS matlab interface
% -----------------------
% author: Thorsten Liebig

info = hdf5info(file);

for n=1:numel(info.GroupHierarchy.Groups)
    if strcmp(info.GroupHierarchy.Groups(n).Name,'/FieldData')
        for m=1:numel(info.GroupHierarchy.Groups(n).Datasets)
            names{m} = info.GroupHierarchy.Groups(n).Datasets(m).Name;
        end
    end
end

hdf_fielddata.names = names;
for n=1:numel(names)
    hdf_fielddata.values{n} = hdf5read(file,names{n});
end