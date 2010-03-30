function invoke_openEMS( opts )

if nargin < 1
    opts = '';
end
% openEMS_opts = [openEMS_opts ' --disable-dumps'];
% openEMS_opts = [openEMS_opts ' --debug-material'];

filename = mfilename('fullpath');
dir = fileparts( filename );
openEMS_Path = [dir '/../../'];

command = [openEMS_Path 'openEMS.sh ' opts];
disp(command);
system(command)
