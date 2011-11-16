function [queue] = Add2Queue(queue,func_name, func_args, varargin)
% function [queue] = Add2Queue(queue,func_name, func_args, varargin)
%
% Use this function to add a funtion to the queue.
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

queue.jobs{jobnum}.finished = 0;

queue.jobs{jobnum}.argsfile = [tempname '.mat'];
save(queue.jobs{jobnum}.argsfile,'func_args');

queue.jobs{jobnum}.nargout = nargout(func_name);
queue.jobs{jobnum}.outargsfile = [tempname '.mat'];

queue.jobs{jobnum}.command = [queue.bin queue.bin_options ' "load(''' queue.jobs{jobnum}.argsfile ''');' ...
    queue.DependPath ...
    '[outargs{1:' num2str(queue.jobs{jobnum}.nargout) '}]=' func_name '(func_args{:});' ...
    'save(''-V7'',''' queue.jobs{jobnum}.outargsfile ''',''outargs'');' ...
    'exit;"'];

[queue.jobs{jobnum}.pid, queue.jobs{jobnum}.filenames] = queue_addProcess( queue.jobs{jobnum}.command );
