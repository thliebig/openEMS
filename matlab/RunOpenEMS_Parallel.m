function [stdout, stderr] = RunOpenEMS_Parallel(Sim_Paths, Sim_Files, opts, Settings, varargin)
% function [stdout, stderr] = RunOpenEMS_Parallel(Sim_Paths, Sim_Files, opts, Settings, varargin)
% 
% Run multiple openEMS simulations in parallel, distributed on multiple 
% machines using a ssh host_list! (currently on Linux only)
% 
% This function relies on InitQueue etc.
% 
% input:
%       Sim_Paths:  cell array of paths to simulate by RunOpenEMS
%       Sim_Files:  filename or cell array of filenames to simulate
%       opts:       openEMS options. see also RunOpenEMS
%       Settings:   use the settings to define multiple host for simulation
%                   e.g.: Settings.SSH.bin ='<path_to_openEMS>/openEMS.sh';
%                         Settings.SSH.host_list = {'list','of','hosts'};
% 
% Note: If no SSH host_list is defined, this function will skip the
%       parallel run and switch back to a default RunOpenEMS!
% 
% See also RunOpenEMS, FindFreeSSH, InitQueue
% 
% openEMS matlab interface
% -----------------------
% author: Thorsten Liebig 2011

pause_queue = 5; %pause between consecutive runs (needed for FindFreeSSH)

skip_parallel = 0;

% currently only supporting linux, run conventional RunOpenEMS
if ~isunix
    warning 'your OS is not supported (Unix only), running default RunOpenEMS';
    skip_parallel = 1;    
end

% in case only one path is given, run conventional RunOpenEMS
if ischar(Sim_Paths)
    warning 'only a single path given, running default RunOpenEMS'
    skip_parallel = 1;
end

% in case SSH.host_list is not defined, run conventional RunOpenEMS
if ~isfield(Settings,'SSH')
    warning 'SSH options missing, running default RunOpenEMS'
    skip_parallel = 1;
elseif ~isfield(Settings.SSH,'host_list')
    warning 'SSH.host_list option missing, running default RunOpenEMS'
    skip_parallel = 1;    
end

if (skip_parallel)
    for n=1:numel(Sim_Paths)
        if iscell(Sim_Files)
            Sim_File = Sim_Files{n};
        else
            Sim_File = Sim_Files;
        end
        RunOpenEMS(Sim_Paths{n}, Sim_Files, opts, Settings)
    end
    stdout = [];
    stderr = [];
    return
end
    
if ~iscell(Sim_Paths)
    error('RunOpenEMS_Parallel:needs a cell array of Sim_Paths to simulate');
end

% get the path to this file
[dir] = fileparts( mfilename('fullpath') );

queue = InitQueue('DependPath',{dir}, varargin{:});

% spawn multiple simulations
for n=1:numel(Sim_Paths)
    if iscell(Sim_Files)
        Sim_File = Sim_Files{n};
    else
        Sim_File = Sim_Files;
    end
    
    queue = Add2Queue(queue,'RunOpenEMS',{Sim_Paths{n}, Sim_File, opts, Settings});
    disp(['openEMS simulation #' int2str(n) ' in directory: ' Sim_Paths{n} ' started!']);
    pause(pause_queue);
end

[queue] = FinishQueue(queue);
   
for n=1:numel(Sim_Paths)
    stdout{n} = queue.jobs{n}.stdout;    
    stderr{n} = queue.jobs{n}.stderr;
end
