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
    openEMS_bin = searchBinary('openEMS.sh', ...
    {[dir filesep '..' filesep '..' filesep], ...  % try devel path
     [dir filesep '..' filesep '..' filesep '..' filesep '..' filesep 'bin' filesep]}); % try (default) install path
else % assume windows
    openEMS_bin = searchBinary('openEMS.exe', [dir filesep '..' filesep '..' filesep]);
end

command = [openEMS_bin ' ' opts];

if ~silent
    if (isunix && nargin>1)
        command = [command ' 2>&1 | tee ' logfile];
    end
else
    command = [command ' > ' logfile ' 2>&1'];
end


exitcode = system(command);
if (exitcode~=0);
	error(['openEMS binary exited with error-code ' num2str(exitcode)]);
end
