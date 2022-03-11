function [queue] = Add2Queue(queue,func_name, func_args, varargin)
% function [queue] = Add2Queue(queue,func_name, func_args, varargin)
%
% Use this function to add a function to the queue.
%
% For more details see: InitQueue
%
% See also: InitQueue, FinishQueue, ResultsQueue, RunOpenEMS
%
% openEMS matlab interface
% -----------------------
% author: Thorsten Liebig

if isfield(queue,'jobs')
    jobnum = numel(queue.jobs)+1;
else
    jobnum = 1;
end

running = numel(queue.jobs_finished) - sum(queue.jobs_finished);

while (running>=queue.maxThreads)
    [queue running] = CheckQueue(queue);
end


if (queue.verbose>=1)
    disp(['Add2Queue: Job #' num2str(jobnum) ' starting...']);
end

queue.jobs_finished(jobnum) = 0;

queue.jobs{jobnum}.argsfile = [tempname '.mat'];
save(queue.jobs{jobnum}.argsfile,'func_args');

queue.jobs{jobnum}.nargout = nargout(func_name);
queue.jobs{jobnum}.outargsfile = [tempname '.mat'];

queue.jobs{jobnum}.command = [queue.bin queue.bin_options ' "load(''' queue.jobs{jobnum}.argsfile ''');' ...
    queue.DependPath ...
    'err=[];' ...
    'try;' ...
    '[outargs{1:' num2str(queue.jobs{jobnum}.nargout) '}]=' func_name '(func_args{:});' ...
    'catch err;outargs=0;end;' ...
    'save(''-V7'',''' queue.jobs{jobnum}.outargsfile ''',''outargs'',''err'');' ...
    'exit;"'];

[queue.jobs{jobnum}.pid, queue.jobs{jobnum}.filenames] = queue_addProcess( queue.jobs{jobnum}.command );
