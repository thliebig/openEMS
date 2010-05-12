function [stdout,stderr] = queue_delProcess( pid, filenames )
% [stdout,stderr] = queue_delProcess( pid, filenames )
%
% if pid == 0, do not kill a process, but clean up files
%
% Sebastian Held <sebastian.held@uni-due.de>
% 12.5.2010

if ~isunix
    error 'your OS is not supported (Unix only)'
end

if pid ~= 0
    alive = queue_checkProcess( pid );

    if alive
        unix( ['kill ' num2str(pid)] );
        alive = queue_checkProcess( pid );
    end
    if alive
        pause(1)
        unix( ['kill ' num2str(pid)] );
        alive = queue_checkProcess( pid );
    end
    if alive
        unix( ['kill -KILL ' num2str(pid)] );
    end
end

if nargin > 1
    if nargout == 1
        [~,stdout] = queue_checkProcess( pid, filenames );
    end
    if nargout == 2
        [~,stdout,stderr] = queue_checkProcess( pid, filenames );
    end

    delete( filenames.stdout );
    delete( filenames.stderr );
end
