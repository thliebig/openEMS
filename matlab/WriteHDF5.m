function WriteHDF5(filename,hdf_fielddata,hdf_mesh)
% function WriteHDF5(filename,hdf_fielddata,hdf_mesh)
%
% input:
%   hdf_fielddata.time
%   hdf_fielddata.names
%   hdf_fielddata.values
%   hdf_mesh.type
%   hdf_mesh.names
%   hdf_mesh.lines
%
% openEMS matlab interface
% -----------------------
% (C) 2010 Sebastian Held <sebastian.held@uni-due.de>
% See also ReadHDF5FieldData ReadHDF5Mesh

isOctave = exist('OCTAVE_VERSION','builtin') ~= 0;
if isOctave
    WriteHDF5_octave(filename,hdf_fielddata,hdf_mesh);
    return
end

writemode = 'overwrite';
if isfield( hdf_fielddata, 'TD' )
    % this is a time domain data set
    time = hdf_fielddata.TD.time;
    for n=1:numel(time)
        name = ['/FieldData/TD/' int2str(n)];
        [details.Location, details.Name] = fileparts(name);
        attribute_details.AttachedTo = name;
        attribute_details.AttachType = 'dataset';
        attribute_details.Name = 'time';
        hdf5write( filename, details, hdf_fielddata.TD.values{n}, ...
                   attribute_details, time(n), ...
                   'WriteMode', writemode );
        writemode = 'append';
    end
end
if isfield( hdf_fielddata, 'FD' )
    % this is a frequency domain data set
    freq = hdf_fielddata.FD.frequency;
    for n=1:numel(freq)
        name = ['/FieldData/FD/f' int2str(n-1) '_real'];
        [details.Location, details.Name] = fileparts(name);
        hdf5write( filename, details, real(hdf_fielddata.FD.values{n}), ...
                   'WriteMode', writemode );
        name = ['/FieldData/FD/f' int2str(n-1) '_imag'];
        [details.Location, details.Name] = fileparts(name);
        hdf5write( filename, details, imag(hdf_fielddata.FD.values{n}), ...
                   'WriteMode', 'append' );
        writemode = 'append';
    end
    name = '/FieldData/FD';
    [details.Location, details.Name] = fileparts(name);
    attribute_details.AttachedTo = name;
    attribute_details.AttachType = 'group';
    attribute_details.Name = 'frequency';
    hdf5write( filename, attribute_details, freq, ...
                   'WriteMode', 'append' );
end

names = hdf_mesh.names; % names is a cell array
for n=1:numel(names)
    [details.Location, details.Name, ext] = fileparts(names{n});
    details.Name = [details.Name ext];
    hdf5write( filename, details, hdf_mesh.lines{n}, ...
               'WriteMode', 'append' );
end



function WriteHDF5_octave(filename,hdf_fielddata,hdf_mesh)
error 'not yet implemented'
