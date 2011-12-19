function [queue] = InitQueue(varargin)
% function [queue] = InitQueue(varargin)
%
% Use this function to initialize a queue to run one or more matlab scripts
% in parallel.
% This can be used to efficiently run an openEMS parameter sweep in parallel
% on multiple remote machines.
%
% Options:
%   DependPath: Add multiple paths, your script may depend on
%   UseOctave:  Enable/Disable octave usage
%   MaxThreads: max. number of parallel executions
%
% Note:
%   - Currently only Linux/Unix is supported
%   - By default Octave is used to spawn parallel functions (saves
%       licenses), but this can be changed by:
%       [queue] = InitQueue('UseOctave', 0);
%       You may need to change this, if your script is not octave compatible
%   - To efficiently run openEMS in parallel, you need to run it on several
%   machines using a SSH.host_list setting --> See also RunOpenEMS
%
% Example:
%   %serial version:
%   for n=1:10
%       % manipulate parameter etc.
%       [result1(n) result2(n)] = Parallel_Func_Name(param1, param2);
%   end
%
%   %parallel version:
%   queue = InitQueue('DependPath',{'/opt/openEMS/CSXCAD/matlab', ...
%                                   '/opt/openEMS/openEMS/matlab'});
%   for n=1:10
%       % manipulate parameter etc.
%       queue = Add2Queue(queue, 'Parallel_Func_Name', {param1, param2});
%   end
%
%   % wait for all to finish
%   [queue] = FinishQueue(queue);
%
%   % retrieve result
%   for n=1:numel(stub_sweep)
%     [result1(n) result2(n)] = ResultsQueue(queue,n);
%   end
%
% See also: Add2Queue, FinishQueue, ResultsQueue, RunOpenEMS,
%           RunOpenEMS_Parallel, FindFreeSSH
%
% openEMS matlab interface
% -----------------------
% author: Thorsten Liebig

if ~isunix
    error 'your OS is not supported (Unix only)'
end

queue.use_octave = exist('OCTAVE_VERSION','builtin') ~= 0;

queue.verbose = 1;

queue.maxThreads = Inf;

% add current path
queue.DependPath = ['addpath(''' pwd ''');'];

for n=1:2:nargin
    if strcmp(varargin{n},'DependPath');
        for m=1:numel(varargin{n+1})
            queue.DependPath = [queue.DependPath 'addpath(''' varargin{n+1}{m} ''');'];
        end
    end
    if strcmp(varargin{n},'UseOctave');
        queue.use_octave = varargin{n+1};
    end
    if strcmp(varargin{n},'MaxThreads');
        queue.maxThreads = varargin{n+1};
    end
end


% set binaries and options
if (queue.use_octave)
    queue.bin = ['export LD_LIBRARY_PATH=""; octave'];
    queue.bin_options = [' --silent --eval'];
else
    queue.bin = [matlabroot '/bin/matlab'];
    queue.bin_options = [' -nodesktop -nosplash  -r'];
end

queue.jobs_finished = [];
