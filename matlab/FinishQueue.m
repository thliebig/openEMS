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

if (queue.verbose>=1)
    disp(['FinishQueue: Waiting for ' num2str(sum(~queue.jobs_finished)) ' of ' num2str(numJobs) ' jobs to finish...']);
end

running = numel(queue.jobs_finished) - sum(queue.jobs_finished);

while sum(running)>0
    [queue running] = CheckQueue(queue, query_time);
end

if (queue.verbose>=1)
    disp(['FinishQueue: All jobs done!'])
end
