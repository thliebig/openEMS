function hdf_fielddata = ReadHDF5FieldData(file)
% function hdf_fielddata = ReadHDF5FieldData(file)
%
% returns:
% % time domain data (if exist)
% hdf_fielddata.TD.time
% hdf_fielddata.TD.names
% hdf_fielddata.TD.values
% hdf_fielddata.TD.DataType (0 --> real value data)
%
% % frequency domain data (if exist)
% hdf_fielddata.FD.frequency
% hdf_fielddata.FD.values
% hdf_fielddata.FD.DataType (0 / 1 --> real / complex value data)
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

if isOctave
    hdf_fielddata = ReadHDF5FieldData_octave(file);
    return
end

info = hdf5info(file);
TD.names = {};
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
        end
        
    end
end

if (numel(TD.names)>0)
    hdf_fielddata.TD=TD;
    hdf_fielddata.TD.DataType = 0; %real value data
    for n=1:numel(hdf_fielddata.TD.names)
        hdf_fielddata.TD.values{n} = double(hdf5read(file,hdf_fielddata.TD.names{n}));
    end
end

% extract FD data
try
    hdf_fielddata.FD.frequency = ReadHDF5Attribute(file,'/FieldData/FD','frequency');
catch err
%     disp(err)
    return
end

for n=1:numel(hdf_fielddata.FD.frequency)
    try
        hdf_fielddata.FD.values{n} = double(hdf5read(file,['/FieldData/FD/f' int2str(n-1) '_real']) + 1i*hdf5read(file,['/FieldData/FD/f' int2str(n-1) '_imag']));
        hdf_fielddata.FD.DataType = 1; %complex value data
    catch
        try
            hdf_fielddata.FD.values{n} = double(hdf5read(file,['/FieldData/FD/f' int2str(n-1)]));
            hdf_fielddata.FD.DataType = 0; %real value data
        catch
            error('openEMS:ReadHDF5FieldData','FD data invalid...')
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
        hdf_fielddata.TD.time(n) = ReadHDF5Attribute(file, hdf_fielddata.TD.names{n},'time');
    end
    hdf_fielddata.TD.DataType = 0; %real value data
end
if isfield(hdf.FieldData,'FD')
    %read FD data
    hdf_fielddata.FD.frequency = ReadHDF5Attribute(file,'/FieldData/FD/','frequency');
    try %try reading complex data
        for n=1:numel(hdf_fielddata.FD.frequency)
            hdf_fielddata.FD.values{n} = double(hdf.FieldData.FD.(['f' int2str(n-1) '_real']) +1i*hdf.FieldData.FD.(['f' int2str(n-1) '_imag']) );
        end
        hdf_fielddata.FD.DataType = 1; %complex value data
    catch
         try %try reading real value data
            for n=1:numel(hdf_fielddata.FD.frequency)
                hdf_fielddata.FD.values{n} = double(hdf.FieldData.FD.(['f' int2str(n-1)]));
            end
            hdf_fielddata.FD.DataType = 0; %real value data
         catch
             error('openEMS:ReadHDF5FieldData','FD data invalid...')
         end
    end
end
