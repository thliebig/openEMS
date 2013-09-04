function invoke_openEMS( opts, logfile, silent )
% function invoke_openEMS( opts, logfile, silent )
%
% internal method to invoke openEMS, use RunOpenEMS instead
%
% See also RunOpenEMS
%
% openEMS matlab interface
% -----------------------
% author: Sebastian Held

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

% opts = [opts ' --disable-dumps'];
% opts = [opts ' --debug-material'];
% opts = [opts ' --debug-boxes'];
% opts = [opts ' --engine=sse'];
% opts = [opts ' --engine=multithreaded'];

filename = mfilename('fullpath');
dir = fileparts( filename );

if isunix
    % <openEMS-path> could be /usr or ~/opt/openEMS etc.
    % assume this file to be in  '<openEMS-path>/share/openEMS/matlab/private/'
    % assume openEMS binary to be in '<openEMS-path>/bin'
    openEMS_Path = [dir filesep '../../../../bin' filesep];
    openEMS_Path = [openEMS_Path 'openEMS.sh'];
else
    openEMS_Path = [dir filesep '../..' filesep];
    openEMS_Path = [openEMS_Path 'openEMS'];
end

command = [openEMS_Path ' ' opts];

if ~silent
    if (isunix && nargin>1)
        command = [command ' 2>&1 | tee ' logfile];
    end
else
    command = [command ' > ' logfile ' 2>&1'];
end

% if ~silent
%     disp( ['invoking openEMS simulator: ' command] );
% end
system(command);
