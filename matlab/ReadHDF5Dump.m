function [field mesh] = ReadHDF5Dump(file, varargin)
%[field mesh] = ReadHDF5Dump(file, varargin)
%
%   Read a hdf5 field dump, including an interpolation and frequency domain
%   transformation.
%
%   For more information about the output, refer to the help of
%   ReadHDF5Mesh and ReadHDF5FieldData
%
%   possible arguments:
%       'Range'             see GetField_Range
%       'Interpolation'     see GetField_Interpolation
%       'SubSampling'       see GetField_SubSampling
%       'Frequency'         see GetField_TD2FD
%       'CloseAlpha': 0 (default) / 1
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
% See also ReadHDF5Mesh ReadHDF5FieldData GetField_Interpolation GetField_SubSampling
% GetField_TD2FD GetField_Range

field = ReadHDF5FieldData(file);
mesh = ReadHDF5Mesh(file);

if (nargin<2)
    return
end

% evaluate arguments in a specific order
for n=1:2:(nargin-1)
    if (strcmp(varargin{n},'Range')==1);
        [field mesh] = GetField_Range(field, mesh, varargin{n+1});
    end
end

for n=1:2:(nargin-1)
    if (strcmp(varargin{n},'SubSampling')==1);
        [field mesh] = GetField_SubSampling(field,mesh,varargin{n+1});
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

for n=1:2:(nargin-1)
    if (strcmp(varargin{n},'CloseAlpha')==1);
        if ((varargin{n+1}==1) && (mesh.type==1) && (range(mesh.lines{2})<2*pi))
            mesh.lines{2}(end+1)=mesh.lines{2}(1)+2*pi;
            if (isfield(field,'TD'))
                for n = 1:numel(field.TD.values)
                    field.TD.values{n}(:,end+1,:,:) =  field.TD.values{n}(:,1,:,:);
                end
            end
            if (isfield(field,'FD'))
                for n = 1:numel(field.FD.values)
                    field.FD.values{n}(:,end+1,:,:) =  field.FD.values{n}(:,1,:,:);
                end
            end
            if (isfield(mesh,'original_indices'))
                if (~isempty(mesh.original_indices))
                    mesh.original_indices{2} = [mesh.original_indices{2} 1];
                end
            else
                mesh.original_indices = {1:numel(mesh.lines{1}),[1:numel(mesh.lines{2}) 1],[1:numel(mesh.lines{3})]};
            end
        end
    end
end

end

function rng = range(x)
    rng = max(x)-min(x);
end