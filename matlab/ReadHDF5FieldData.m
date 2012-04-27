function hdf_fielddata = ReadHDF5FieldData(file)
% function hdf_fielddata = ReadHDF5FieldData(file)
%
% returns:
% % time domain data (if exist)
% hdf_fielddata.TD.time
% hdf_fielddata.TD.names
% hdf_fielddata.TD.values
%
% % frequency domain data (if exist)
% hdf_fielddata.FD.time
% hdf_fielddata.FD.names
% hdf_fielddata.FD.values
%
% example: values of timestep 12:
% hdf_fielddata.TD.values{12}: array (x,y,z,polarization)
%
% plot z-field component along y-direction for timestep 12:
% plot( hdf_fielddata.TD.values{12}(1,:,1,3) )
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
                    for a = 1:numel(info.GroupHierarchy.Groups(n).Groups(nGroup).Datasets(m).Attributes)
                        str = regexp(info.GroupHierarchy.Groups(n).Groups(nGroup).Datasets(m).Attributes(a).Name,'\w/*\w*','match');
                        TD.(str{end})(m) = double(info.GroupHierarchy.Groups(n).Groups(nGroup).Datasets(m).Attributes(a).Value);
                    end
                end
            end
            %search and read FD data
            if strcmp(info.GroupHierarchy.Groups(n).Groups(nGroup).Name,'/FieldData/FD')
                for m=1:numel(info.GroupHierarchy.Groups(n).Groups(nGroup).Datasets)
                    FD.names{m} = info.GroupHierarchy.Groups(n).Groups(nGroup).Datasets(m).Name;
                    for a = 1:numel(info.GroupHierarchy.Groups(n).Groups(nGroup).Datasets(m).Attributes)
                        str = regexp(info.GroupHierarchy.Groups(n).Groups(nGroup).Datasets(m).Attributes(a).Name,'\w/*\w*','match');
                        FD.(str{end})(m) = double(info.GroupHierarchy.Groups(n).Groups(nGroup).Datasets(m).Attributes(a).Value);
                    end
                end
            end
        end
        
    end
end

if (numel(TD.names)>0)
    hdf_fielddata.TD=TD;
    for n=1:numel(hdf_fielddata.TD.names)
        hdf_fielddata.TD.values{n} = double(hdf5read(file,hdf_fielddata.TD.names{n}));
    end
end

if (numel(FD.names)>0)
    hdf_fielddata.FD=FD;
    Nr_freq = numel(FD.names);
    for n=1:Nr_freq
        name = ['/FieldData/FD/f' int2str(n-1) '_real'];
        ind = find(strcmp(FD.names,name));
        if isempty(ind)
            ind = find(strcmp(FD.names,['/FieldData/FD/f' int2str(n-1)]));
            if ~isempty(ind)
                hdf_fielddata.FD.values{n} = double(hdf5read(file,FD.names{ind}));
            end
        else
            hdf_fielddata.FD.values{n} = double(hdf5read(file,FD.names{ind}));
            name = ['/FieldData/FD/f' int2str(n-1) '_imag'];
            ind = find(strcmp(FD.names,name));
            hdf_fielddata.FD.values{n} = hdf_fielddata.FD.values{n} + 1j*double(hdf5read(file,FD.names{ind}));
        end
    end
end

function hdf_fielddata = ReadHDF5FieldData_octave(file)
hdf = load( '-hdf5', file );
if ~isfield(hdf,'FieldData')
    error('no field data found')
end
if isfield(hdf.FieldData,'TD')
    %read TD data
    hdf_fielddata_names = fieldnames(hdf.FieldData.TD);
    for n=1:numel(hdf_fielddata_names)
        hdf_fielddata.TD.values{n} = hdf.FieldData.TD.(hdf_fielddata_names{n});
        hdf_fielddata.TD.names{n} = ['/FieldData/TD/' hdf_fielddata_names{n}(2:end)];
        hdf_fielddata.TD.time(n) = h5readatt_octave(file, hdf_fielddata.TD.names{n},'time');
    end
end
if isfield(hdf.FieldData,'FD')
    %read FD data
    hdf_fielddata_names = fieldnames(hdf.FieldData.FD);
    for n=1:numel(hdf_fielddata_names)/2
        hdf_fielddata.FD.values{n} = hdf.FieldData.FD.(hdf_fielddata_names{2*n}) + 1j*hdf.FieldData.FD.(hdf_fielddata_names{2*n-1});
        hdf_fielddata.FD.names{n} = ['/FieldData/FD/' hdf_fielddata_names{2*n-1}(1:end-5)];
        hdf_fielddata.FD.frequencies(n) = h5readatt_octave(file,['/FieldData/FD/' hdf_fielddata_names{2*n}],'frequency');
    end
end
