function path = openEMS_resource_path(varargin)
% path = openEMS_resource_path(varargin)
%
% Absolute path of a data file shipped with openEMS, e.g.
%
%   phantom = openEMS_resource_path('phantoms', 'phantom_head_298MHz.h5');
%
% The path is derived from the location of this function, not from the
% caller, so scripts using it work from any directory.
%
% arguments:
%   varargin : path components below the openEMS resources folder
%
% openEMS matlab interface
% -----------------------
% author: Thorsten Liebig

% <interface>/../resources: holds in the source tree, in an installation
% (share/openEMS/) and in the windows package alike.
root = fullfile(fileparts(mfilename('fullpath')), '..', 'resources');
path = fullfile(root, varargin{:});

if ~exist(path, 'file')
    error('openEMS:resource_not_found', ...
          'openEMS resource not found: %s', path);
end
