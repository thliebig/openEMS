function [field mesh] = ReadHDF5Dump(file, varargin)
%[field mesh] = ReadHDF5Dump(file, varargin)
% 
%   Read a hdf5 field dump, including an interpolation and frequency domain
%   transformation.
% 
%   possible arguments:
%       'Range'             see GetField_Range
%       'Interpolation'     see GetField_Interpolation
%       'Frequency'         see GetField_TD2FD
%   
%   example:
%   [field mesh] = ReadHDF5Dump('Et.h5');
%   or
%   [field mesh] = ReadHDF5Dump('Et.h5','Range',{[0 100],[-20 20],[50 90]});
%   or
%   [field mesh] = ReadHDF5Dump('Et.h5','Interpolation',[21 1 101],'Frequency',300e6);
%
% openEMS matlab interface
% -----------------------
% author: Thorsten Liebig
%
% See also ReadHDF5Mesh ReadHDF5FieldData GetField_Interpolation
% GetField_TD2FD GetField_Range

field = ReadHDF5FieldData(file);
mesh = ReadHDF5Mesh(file);

if (nargin>1)
    for n=1:2:(nargin-1)
        if (strcmp(varargin{n},'Range')==1);
            [field mesh] = GetField_Range(field, mesh, varargin{n+1});
        end
    end
    
    for n=1:2:(nargin-1)
        if (strcmp(varargin{n},'Interpolation')==1);
            [field mesh] = GetField_Interpolation(field,mesh,varargin{n+1});
        end
    end
    
    for n=1:2:(nargin-1)
        if (strcmp(varargin{n},'Frequency')==1);
            field = GetField_TD2FD(field,varargin{n+1});
        end
    end
end