function invoke_openEMS( opts, logfile, silent )
% function invoke_openEMS( opts, logfile, silent )
%
% internal method to invoke openEMS, use RunOpenEMS instead
%
% See also RunOpenEMS
%
% openEMS matlab interface
% -----------------------
% author: Sebastian Held, Thorsten Liebig

if nargin < 1
    error 'specify the xml file to simulate'
end
if nargin < 3
    silent = 0;
end
if (nargin < 2) || isempty(logfile)
    if isunix
        logfile = '/dev/null';
    else
        logfile = 'nul:';
    end
end

filename = mfilename('fullpath');
dir = fileparts( filename );

if isunix
    % try development path
    openEMS_bin = [dir filesep '../..' filesep 'openEMS.sh'];
    if (~exist(openEMS_bin,'file'))
        % fallback to install path
        openEMS_bin = [dir filesep '../../../../bin' filesep 'openEMS.sh'];
    end
else
    openEMS_bin = [dir filesep '../..' filesep];
    openEMS_bin = [openEMS_bin 'openEMS'];
end

if (~exist(openEMS_bin,'file'))
    error('openEMS:invoke_openEMS',['Binary not found: ' openEMS_bin]);
end

command = [openEMS_bin ' ' opts];

if ~silent
    if (isunix && nargin>1)
        command = [command ' 2>&1 | tee ' logfile];
    end
else
    command = [command ' > ' logfile ' 2>&1'];
end

system(command);
