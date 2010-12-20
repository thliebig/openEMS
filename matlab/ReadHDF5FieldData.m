function hdf_fielddata = ReadHDF5FieldData(file)
% function hdf_fielddata = ReadHDF5FieldData(file)
%
% returns:
% hdf_fielddata.time
% hdf_fielddata.names
% hdf_fielddata.values
%
% example: values of timestep 12:
% hdf_fielddata.values{12}: array (x,y,z,polarization)
%
% plot z-field component along y-direction for timestep 12:
% plot( hdf_fielddata.values{12}(1,:,1,3) )
%
% openEMS matlab interface
% -----------------------
% author: Thorsten Liebig
%
% See also ReadHDF5Mesh ReadHDF5Dump

isOctave = exist('OCTAVE_VERSION','builtin') ~= 0;
if isOctave
    hdf_fielddata = ReadHDF5FieldData_octave(file);
    return
end

info = hdf5info(file);
TD.names = {};
FD.names = {};
hdf_fielddata = [];

for n=1:numel(info.GroupHierarchy.Groups)
    if strcmp(info.GroupHierarchy.Groups(n).Name,'/FieldData')
        %found /FieldData, look for either TD or FD data
        for nGroup=1:numel(info.GroupHierarchy.Groups(n).Groups)
            %search and read TD data
            if strcmp(info.GroupHierarchy.Groups(n).Groups(nGroup).Name,'/FieldData/TD')
                for m=1:numel(info.GroupHierarchy.Groups(n).Groups(nGroup).Datasets)
                    TD.names{m} = info.GroupHierarchy.Groups(n).Groups(nGroup).Datasets(m).Name;
                    hdf_fielddata.TD.time(m) = double(info.GroupHierarchy.Groups(n).Groups(nGroup).Datasets(m).Attributes.Value);
                end
            end
            %search and read FD data
            if strcmp(info.GroupHierarchy.Groups(n).Groups(nGroup).Name,'/FieldData/FD')
                for m=1:numel(info.GroupHierarchy.Groups(n).Groups(nGroup).Datasets)
                    FD.names{m} = info.GroupHierarchy.Groups(n).Groups(nGroup).Datasets(m).Name;
                    FD.freq(m) = double(info.GroupHierarchy.Groups(n).Groups(nGroup).Datasets(m).Attributes.Value);
                end
            end
        end
        
    end
end

if (numel(TD.names)>0)
    for n=1:numel(TD.names)
        hdf_fielddata.TD.values{n} = double(hdf5read(file,TD.names{n}));
    end
end

if (numel(FD.names)>0)
    Nr_freq = numel(FD.names)/2;
    for n=1:Nr_freq
        name = ['/FieldData/FD/f' int2str(n-1) '_real'];
        ind = find(strcmp(FD.names,name));
        hdf_fielddata.FD.values{n} = double(hdf5read(file,FD.names{ind}));
        hdf_fielddata.FD.freq(n) = FD.freq(ind);
        name = ['/FieldData/FD/f' int2str(n-1) '_imag'];
        ind = find(strcmp(FD.names,name));
        hdf_fielddata.FD.values{n} = hdf_fielddata.FD.values{n} + 1j*double(hdf5read(file,FD.names{ind}));
    end
end


function hdf_fielddata = ReadHDF5FieldData_octave(file)
hdf = load( '-hdf5', file );
hdf_fielddata.names = fieldnames(hdf.FieldData);
for n=1:numel(hdf_fielddata.names)
    hdf_fielddata.time(n) = str2double(hdf_fielddata.names{n}(2:end));
    hdf_fielddata.values{n} = hdf.FieldData.(hdf_fielddata.names{n});
end
