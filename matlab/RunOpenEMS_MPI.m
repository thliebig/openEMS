function RunOpenEMS_MPI(Sim_Path, Sim_File, opts, Settings)
% function RunOpenEMS_MPI(Sim_Path, Sim_File, NrProc,  opts, Settings)
%
% Run an openEMS simulation with MPI support
% 
% % mpi binary path on all nodes needed
% Settings.MPI.Binary = '~/devel/openEMS/openEMS_MPI';
% % number of processes to run
% Settings.MPI.NrProc = 2;
% 
% % optional:
% % define a hostfile and number of host to boot the mpd daemon:
% Settings.MPI.HostFile = '/home/thorsten/ate-pc9x.hosts';
% Settings.MPI.TotalNum = 2;
% 
% RunOpenEMS_MPI(Sim_Path, Sim_File, NrProc, opts, Settings)
%
% See also WriteOpenEMS, RunOpenEMS
%
% openEMS matlab interface
% -----------------------
% author: Thorsten Liebig

if (isunix ~= 1)
    error 'MPI version of openEMS currently only available using Linux'
end

if nargin < 4
    error 'missing arguments: specify the Sim_Path, Sim_file, opts and Settings...'
end

NrProc = Settings.MPI.NrProc;

if (NrProc<2)
    error('openEMS:RunOpenEMS_MPI','MPI number of processes to small...');
end

if ~isfield(Settings,'MPI')
    error('openEMS:RunOpenEMS_MPI','MPI settings not found...');
end

if isfield(Settings.MPI,'HostFile')
    [status, result] = unix(['mpdboot -v --file=' Settings.MPI.HostFile ' --totalnum=' int2str(Settings.MPI.TotalNum)]);
    if (status~=0)
        disp(result);
        error('openEMS:RunOpenEMS','mpdboot failed to boot mpi daemon!');
    end
else
    [status, result] = unix(['mpdboot -v']);
    if (status~=0)
        disp(result);
        error('openEMS:RunOpenEMS','mpdboot failed to boot mpi daemon!');
    end
end

savePath = pwd;
cd(Sim_Path);

scp_options = '-C -o "PasswordAuthentication no" -o "StrictHostKeyChecking no"';
ssh_options = [scp_options ' -x'];

if isfield(Settings.MPI,'HostFile')
    [status, result] = unix(['mpirun -machinefile ' Settings.MPI.HostFile ' -n  ' int2str(NrProc) '  hostname']);
else
    [status, result] = unix(['mpirun -n ' int2str(NrProc) '  hostname']);
end

if (status~=0)
    disp(result);
    error('openEMS:RunOpenEMS',['mpirun failed ...']);
end

Remote_Nodes = regexp(result, '([^ \n][^\n]*)', 'match'); %get the names of all mpi nodes
Remote_Nodes = unique(Remote_Nodes);

for n=1:numel(Remote_Nodes)
    remote_name = Remote_Nodes{n};
    
    if (n==1)
        [status, result] = unix(['ssh ' ssh_options ' ' remote_name ' "mktemp -d /tmp/openEMS_MPI_XXXXXXXXXXXX"']);
        if (status~=0)
            disp(result);
            error('openEMS:RunOpenEMS','mktemp failed to create tmp directory!');
        end
        work_path = strtrim(result); %remove tailing \n
    else
        [status, result] = unix(['ssh ' ssh_options ' ' remote_name ' "mkdir ' work_path '"']);       
        if (status~=0)
            disp(result);
            error('openEMS:RunOpenEMS',['mkdir failed to create tmp directory on remote ' remote_name ' !']);
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

disp(['Running remote openEMS_MPI in working dir: ' work_path]);

if ~isfield(Settings.MPI,'GlobalArgs')
    Settings.MPI.GlobalArgs = '';
end

if isfield(Settings.MPI,'HostFile')
    [status]  = system(['LD_LIBRARY_PATH= mpirun -machinefile ' Settings.MPI.HostFile ' -l -n ' int2str(NrProc) ' -wdir ' work_path ' ' Settings.MPI.Binary ' ' Sim_File ' ' opts ' ' append_unix]);
else
    [status]  = system(['LD_LIBRARY_PATH= mpirun ' Settings.MPI.GlobalArgs ' -n ' int2str(NrProc) ' -wdir ' work_path ' ' Settings.MPI.Binary ' ' Sim_File ' ' opts ' ' append_unix]);
end
if (status~=0)
    error('openEMS:RunOpenEMS','mpirun openEMS failed!');
end

disp( 'Remote simulation done... copying back results and cleaning up...' );

if (strncmp(work_path,'/tmp/',5)~=1) % savety precaution...
    error('openEMS:RunOpenEMS','working path invalid for deletion');
end
    
for n=1:numel(Remote_Nodes)
    remote_name = Remote_Nodes{n};
    disp(['Copy data from remote node: ' remote_name]);
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
