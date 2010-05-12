function [pid,filenames] = queue_addProcess( command )
% [pid,filenames] = queue_addProcess( command )
%
% Sebastian Held <sebastian.held@uni-due.de>
% 12.5.2010

if ~isunix
    error 'your OS is not supported (Unix only)'
end

if nargout > 1
    filenames.stdout = tempname;
    filenames.stderr = tempname;
else
    filenames.stdout = '/dev/null';
    filenames.stderr = '/dev/null';
end

cmd = ['(' command ') >' filenames.stdout ' 2>' filenames.stderr ' & echo $!' ];
[~,result] = unix( cmd );

pid = str2double(result);
