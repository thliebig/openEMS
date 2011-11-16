function [queue] = FinishQueue(queue, query_time)
% function [queue] = FinishQueue(queue, <query_time>)
%
% Wait for the given queue to finish.
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
    return
end

if (nargin<2)
    query_time = 5;
end

numJobs = numel(queue.jobs);

for n=1:numel(numJobs)
    is_done = queue.jobs{n}.finished;
end

if (queue.verbose>=1)
    disp(['FinishQueue: Waiting for ' num2str(sum(~is_done)) ' of ' num2str(numJobs) ' jobs to finish...']);
end

while sum(is_done)<numJobs
    pause(query_time);

    for n=1:numel(numJobs)
        if (is_done(n)==0)
            if (queue_checkProcess( queue.jobs{n}.pid, queue.jobs{n}.filenames)==0)
                queue.jobs{n}.finished=1;
                load(queue.jobs{n}.outargsfile);
                queue.jobs{n}.outargs = outargs;

                % read in output and cleanup
                [queue.jobs{n}.stdout,queue.jobs{n}.stderr] = queue_delProcess( queue.jobs{n}.pid, queue.jobs{n}.filenames );

                % cleanup
                delete( queue.jobs{n}.argsfile );
                clear queue.jobs{n}.argsfile;
                delete( queue.jobs{n}.outargsfile );
                clear queue.jobs{n}.outargsfile;

                is_done(n) = 1;

                if (queue.verbose>=1)
                    disp(['FinishQueue: Job #' num2str(n) ' is finished!']);
                end
            end
        end
    end
end

if (queue.verbose>=1)
    disp(['FinishQueue: All jobs done!'])
end
