function RunOpenEMS_MPI(Sim_Path, Sim_File, opts, Settings)
% function RunOpenEMS_MPI(Sim_Path, Sim_File, NrProc,  opts, Settings)
%
% Run an openEMS simulation with MPI support
%
% % mpi binary path on all nodes needed
% Settings.MPI.Binary = '/opt/openEMS/openEMS';
% % number of processes to run
% Settings.MPI.NrProc = 3;
% % define the mpi hosts :
% Settings.MPI.Hosts = {'host1','host2','host3'};
%
% RunOpenEMS(Sim_Path, Sim_File, NrProc, opts, Settings)
%
% See also SetupMPI, WriteOpenEMS, RunOpenEMS
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

savePath = pwd;
cd(Sim_Path);

scp_options = '-C -o "PasswordAuthentication no" -o "StrictHostKeyChecking no"';
ssh_options = [scp_options ' -x'];

if isfield(Settings.MPI,'Hosts')
    Remote_Nodes = Settings.MPI.Hosts;
    HostList = '';
    for n=1:numel(Remote_Nodes)
        remote_name = Remote_Nodes{n};
        
        if (n==1)
            [status, result] = unix(['ssh ' ssh_options ' ' remote_name ' "mktemp -d /tmp/openEMS_MPI_XXXXXXXXXXXX"']);
            if (status~=0)
                disp(result);
                error('openEMS:RunOpenEMS','mktemp failed to create tmp directory!');
            end
            work_path = strtrim(result); %remove tailing \n
            HostList = remote_name;
        else
            [status, result] = unix(['ssh ' ssh_options ' ' remote_name ' "mkdir ' work_path '"']);
            if (status~=0)
                disp(result);
                error('openEMS:RunOpenEMS',['mkdir failed to create tmp directory on remote ' remote_name ' !']);
            end
            HostList = [HostList ',' remote_name];
        end
        
        [stat, res] = unix(['scp ' scp_options ' * ' remote_name ':' work_path '/']);
        if (stat~=0)
            disp(res);
            error('openEMS:RunOpenEMS',['scp to remote ' remote_name ' failed!']);
        end
    end
end


%run openEMS (with log file if requested)
if isfield(Settings,'LogFile')
    append_unix = [' 2>&1 | tee ' Settings.LogFile];
else
    append_unix = [];
end

if ~isfield(Settings.MPI,'GlobalArgs')
    Settings.MPI.GlobalArgs = '';
end

if isfield(Settings.MPI,'Hosts')
    disp(['Running remote openEMS_MPI in working dir: ' work_path]);
    [status]  = system(['mpiexec ' Settings.MPI.GlobalArgs ' -host ' HostList ' -n ' int2str(NrProc) ' -wdir ' work_path ' ' Settings.MPI.Binary ' ' Sim_File ' ' opts ' ' append_unix]);
else
    disp('Running local openEMS_MPI');
    [status]  = system(['mpiexec ' Settings.MPI.GlobalArgs ' -n ' int2str(NrProc) ' ' Settings.MPI.Binary ' ' Sim_File ' ' opts ' ' append_unix]);
end

if (status~=0)
    error('openEMS:RunOpenEMS','mpirun openEMS failed!');
end

if isfield(Settings.MPI,'Hosts')
    disp( 'Remote simulation done... copying back results and cleaning up...' );
    
    if (strncmp(work_path,'/tmp/',5)~=1) % safety precaution...
        error('openEMS:RunOpenEMS','working path invalid for deletion');
    end
    
    for n=1:numel(Remote_Nodes)
        remote_name = Remote_Nodes{n};
        disp(['Copy data from remote node: ' remote_name]);
        [stat, res] = unix(['scp -r ' scp_options ' ' remote_name ':' work_path '/* ''' pwd '''/']);
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
end

cd(savePath);
