function invoke_openEMS( opts , logfile)
% function invoke_openEMS( opts )
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

% opts = [opts ' --disable-dumps'];
% opts = [opts ' --debug-material'];
% opts = [opts ' --debug-boxes'];
% opts = [opts ' --engine=sse'];
% opts = [opts ' --engine=multithreaded'];

filename = mfilename('fullpath');
dir = fileparts( filename );
openEMS_Path = [dir filesep '..' filesep];
    
if isunix
	openEMS_Path = [openEMS_Path 'openEMS.sh'];
else
	openEMS_Path = [openEMS_Path 'openEMS'];
end

command = [openEMS_Path ' ' opts];

if (isunix && nargin>1)
    command = [command ' 2>&1 | tee ' logfile];
end

disp( ['invoking openEMS simulator: ' command] );
system(command);
