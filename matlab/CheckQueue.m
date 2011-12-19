function [queue running] = CheckQueue(queue, query_time)
% function [queue running] = CheckQueue(queue, <query_time>)
%
% Check the given queue for finished tasks.
%
% Parameter:
%   query_time (optional): time interval to check for finished tasks
%                          (in seconds, default is 5)
%
% For more details see: InitQueue
%
% See also: InitQueue, ResultsQueue, Add2Queue, RunOpenEMS
%
% openEMS matlab interface
% -----------------------
% author: Thorsten Liebig

if ~isfield(queue,'jobs')
    running = 0;
    return
end

if (nargin<2)
    query_time = 5;
end

numJobs = numel(queue.jobs);

pause(query_time);

for n=1:numJobs
    if (queue.jobs_finished(n)==0)
        if (queue_checkProcess( queue.jobs{n}.pid, queue.jobs{n}.filenames)==0)
            queue.jobs_finished(n)=1;
            load(queue.jobs{n}.outargsfile);
            if ~isempty(err)
                disp(['Job with number ' num2str(n) ' failed to execute: Error message:']);
                error(['CheckQueue:' err.message]);
            end
            queue.jobs{n}.outargs = outargs;

            % read in output and cleanup
            [queue.jobs{n}.stdout,queue.jobs{n}.stderr] = queue_delProcess( queue.jobs{n}.pid, queue.jobs{n}.filenames );

            % cleanup
            delete( queue.jobs{n}.argsfile );
            clear queue.jobs{n}.argsfile;
            delete( queue.jobs{n}.outargsfile );
            clear queue.jobs{n}.outargsfile;

            queue.jobs_finished(n) = 1;

            if (queue.verbose>=1)
                disp(['CheckQueue: Job #' num2str(n) ' is finished!']);
            end
        end
    end
end

running = numel(queue.jobs_finished) - sum(queue.jobs_finished);
