function [alive,stdout,stderr] = queue_checkProcess( pid, filenames )
% [alive,stdout,stderr] = queue_checkProcess( pid )
%
% Sebastian Held <sebastian.held@uni-due.de>
% 12.5.2010

if ~isunix
    error 'your OS is not supported (Unix only)'
end

if nargout > 1
    fid = fopen( filenames.stdout );
    stdout = fread(fid, '*char')';
    fclose(fid);
end
if nargout > 2
    fid = fopen( filenames.stderr );
    stderr = fread(fid, '*char')';
    fclose(fid);
end

cmd = ['ps --no-headers -p' num2str(pid) ];
[status,~] = unix( cmd );

alive = (status == 0);
