function [varargout] = ResultsQueue(queue, n)
% function [varargout] = ResultsQueue(queue, n)
%
% Use this function to retrieve the results from a finished queue.
%
% For more details see: InitQueue
%
% See also: InitQueue, FinishQueue, Add2Queue, RunOpenEMS
%
% openEMS matlab interface
% -----------------------
% author: Thorsten Liebig

if n>numel(queue.jobs)
    error 'ResultsQueue:job is missing'
end

if (nargout>numel(queue.jobs{n}.outargs))
    error 'not enough job output arguments'
end

for k=1:numel(queue.jobs{n}.outargs)
    varargout{k} = queue.jobs{n}.outargs{k};
end
