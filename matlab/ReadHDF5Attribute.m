function attr = ReadHDF5Attribute(file, groupname, attr_name)
% attr = ReadHDF5Attribute(file, groupname, attr_name)
%
% internal function for openEMS to read hdf5 attributes
%
% See also: ReadHDF5ComplexData
%
% openEMS Matlab/Octave interface
% -----------------------
% author: Thorsten Liebig, 2012


if isOctave
    if (exist('h5readatt_octave')==0)
        warning('openEMS:ReadHDF5Attribute','function "h5readatt_octave" not found, trying to run "setup"');
        try
            setup
        catch
            error('openEMS:ReadHDF5Attribute','running "setup" failed...');
        end
    end
    attr = double(h5readatt_octave(file,groupname,attr_name));
else
    %check for different matlab versions
    if verLessThan('matlab','7.9')
        attr = double(hdf5read(file,[groupname '/' attr_name]));
    elseif verLessThan('matlab','7.12')
        attr = double(hdf5read(file,groupname,attr_name));
    else
        attr = double(h5readatt(file,groupname,attr_name));
    end
    
end