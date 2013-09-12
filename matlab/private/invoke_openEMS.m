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
    openEMS_bin = searchBinary('openEMS.sh',[dir filesep '..' filesep '..' filesep]);
else % assume windows
    openEMS_bin = searchBinary('openEMS.exe',[dir filesep '..' filesep '..' filesep]);
end

if (isempty(openEMS_bin))
    error('openEMS:invoke_openEMS', 'openEMS binary not found!');
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
