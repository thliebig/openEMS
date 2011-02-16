function RunOpenEMS_MPI(Sim_Path, Sim_File, NrProc, opts, Settings, copy_bin)
% function RunOpenEMS_MPI(Sim_Path, Sim_File, NrProc,  opts, Settings, copy_bin)
%
% Run an openEMS simulation with MPI support
%
% example:
% Sim_Path = 'MySimPath';
% Sim_File = 'helix.xml'; %should be created by WriteOpenEMS
% 
% NrProc = 2; % set the number of processes to start
% 
% % mpi binary path on all nodes needed
% Settings.MPI_Binary = '~/devel/openEMS/openEMS_MPI';
% 
% opts = '--engine=MPI';
% 
% optional: 
% Settings.LogFile = 'openEMS.log'
% Settings.Silent  = 0
%
% copy_bin = 1; %copy openEMS binary to <Settings.MPI_Binary> on all nodes (mainly for developers)
%
% RunOpenEMS(Sim_Path,Sim_File, NrProc, opts, Settings)
%
% See also WriteOpenEMS, RunOpenEMS
%
% openEMS matlab interface
% -----------------------
% author: Thorsten Liebig

if (isunix ~= 1)
    error 'MPI version of openEMS currently only available using Linux'
end

if nargin < 5
    error 'missing arguments: specify the Sim_Path, Sim_file, Nodes and NrProc...'
end

if (NrProc<2)
    warning('openEMS:RunOpenEMS_MPI','MPI number of processes to small... running non-MPI openEMS');
    RunOpenEMS(Sim_Path,Sim_File,opts,Settings);
    return;
end

if nargin < 6
    copy_bin = 0;
end

savePath = pwd;
cd(Sim_Path);

% setup tmp directory on host machine 

[status, result] = unix('mktemp -d /tmp/openEMS_MPI_XXXXXXXXXXXX');
if (status~=0)
    disp(result);
    error('openEMS:RunOpenEMS','mktemp failed to create tmp directory!');
end

work_path = strtrim(result); %remove tailing \n

disp(['Running remote openEMS_MPI in working dir: ' work_path]);

%copy openEMS & all simulation files to host
if (copy_bin>0)
    filename = mfilename('fullpath');
    dir = fileparts( filename );
    openEMS_Path = [dir filesep '..' filesep];

    [stat, res] = unix(['cp ' openEMS_Path 'openEMS ' Settings.MPI_Binary]);
    if (stat~=0)
        disp(res);
        error('openEMS:RunOpenEMS','host copy failed!');
    end
end

[stat, res] = unix(['cp * ' work_path '/']);
if (stat~=0)
    disp(res);
    error('openEMS:RunOpenEMS','host copy failed!');
end

scp_options = '-C -o "PasswordAuthentication no" -o "StrictHostKeyChecking no"';
ssh_options = [scp_options ' -x'];

[status, result] = unix(['mpirun -n ' int2str(NrProc) '  hostname']);
if (status~=0)
    disp(result);
    error('openEMS:RunOpenEMS',['mpirun failed ...']);
end

[status, LocalNode] = unix('hostname'); %name of local host
Remote_Nodes = regexp(result, '([^ \n][^\n]*)', 'match'); %get the names of all mpi nodes
    
Remote_Nodes = setdiff(Remote_Nodes,LocalNode); %remove local host from node list

for n=1:numel(Remote_Nodes)
    remote_name = Remote_Nodes{n};
    
 	[status, result] = unix(['ssh ' ssh_options ' ' remote_name ' "mkdir ' work_path '"']);
    if (status~=0)
        disp(result);
        error('openEMS:RunOpenEMS',['mkdir failed to create tmp directory on remote ' remote_name ' !']);
    end
    
    %copy openEMS & all simulation files to the ssh host
    if (copy_bin>0)
        [stat, res] = unix(['scp ' scp_options ' ' Settings.MPI_Binary ' ' remote_name ':' Settings.MPI_Binary]);
        if (stat~=0)
            disp(res);
            error('openEMS:RunOpenEMS',['scp to remote ' remote_name ' failed!']);
        end
    end
    
    [stat, res] = unix(['scp ' scp_options ' * ' remote_name ':' work_path '/']);
    if (stat~=0)
        disp(res);
        error('openEMS:RunOpenEMS',['scp to remote ' remote_name ' failed!']);
    end
end


%run openEMS (with log file if requested)
if isfield(Settings,'LogFile')
    append_unix = [' 2>&1 | tee ' Settings.LogFile];
else
    append_unix = [];
end

status = system(['LD_LIBRARY_PATH= mpirun -l -n ' int2str(NrProc) ' -wdir ' work_path ' ' Settings.MPI_Binary ' ' Sim_File ' ' opts ' ' append_unix]);
if (status~=0)
    disp(result);
    error('openEMS:RunOpenEMS','mpirun openEMS failed!');
end

disp( 'Remote simulation done... copying back results and cleaning up...' );

if (strncmp(work_path,'/tmp/',5)~=1) % savety precaution...
    error('openEMS:RunOpenEMS','working path invalid for deletion');
end

%copy back all results
[stat, res] = unix(['cp -r ' work_path '/* ' pwd '/']);
if (stat~=0);
    disp(res);
    error('openEMS:RunOpenEMS','host cp failed!');
end

%cleanup
[stat, res] = unix([' rm -r ' work_path]);
if (stat~=0);
    disp(res);
    warning('openEMS:RunOpenEMS','host cleanup failed!');
end
    
for n=1:numel(Remote_Nodes)
    remote_name = Remote_Nodes{n};
    [stat, res] = unix(['scp -r ' scp_options ' ' remote_name ':' work_path '/* ' pwd '/']);
    if (stat~=0);
        disp(res);
        error('openEMS:RunOpenEMS','remote scp failed!');
    end

    %cleanup
    [stat, res] = unix(['ssh ' ssh_options ' ' remote_name ' rm -r ' work_path]);
    if (stat~=0);
        disp(res);
        warning('openEMS:RunOpenEMS','remote cleanup failed!');
    end
end
    
cd(savePath);
