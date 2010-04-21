function invoke_openEMS( opts )

if nargin < 1
    opts = '';
end
% opts = [opts ' --disable-dumps'];
% opts = [opts ' --debug-material'];
% opts = [opts ' --engine=multithreaded'];
% opts = [opts ' --engine=sse'];

filename = mfilename('fullpath');
dir = fileparts( filename );
openEMS_Path = [dir '/../../'];

command = [openEMS_Path 'openEMS.sh ' opts];
disp(command);
system(command);
