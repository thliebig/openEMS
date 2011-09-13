function RunOpenEMS(Sim_Path, Sim_File, opts, Settings)
% function RunOpenEMS(Sim_Path, Sim_File, opts, Settings)
%
% Run an openEMS simulation
%
% example:
% Sim_Path = 'MySimPath';
% Sim_File = 'helix.xml'; %should be created by WriteOpenEMS
% opts = '--engine=fastest';
% 
% optional:  
% (ssh only on unix with working ssh client)
% Settings.SSH.host = '<hostname or ip>'
% Settings.SSH.bin = '<path_to_openEMS>/openEMS.sh'
% ssh optional: 
% Settings.SSH.host_list = {'list','of','hosts'}; %searches for a free host
%
% optional MPI:
% Settings.MPI.xxx --> help RunOpenEMS_MPI
% 
% Settings.LogFile = 'openEMS.log'
% Settings.Silent  = 0
%
% RunOpenEMS(Sim_Path,Sim_File,opts,Settings)
%
% See also WriteOpenEMS FindFreeSSH
%
% openEMS matlab interface
% -----------------------
% author: Thorsten Liebig

if nargin < 2
    error 'specify the Sim_Path and Sim_file to simulate'
end

if nargin < 3
    opts = '';
end

if (nargin<4)
    Settings = [];
end

if (isfield(Settings,'MPI') && isunix)
    if (Settings.MPI.NrProc>1)
        RunOpenEMS_MPI(Sim_Path, Sim_File, opts, Settings);
        return;
    end
end

savePath = pwd;
cd(Sim_Path);
    
if (isfield(Settings,'SSH') && isunix)
    % ssh options: no X forwarding; no password prompt (use pub keys!); no
    % host checking
    scp_options = '-C -o "PasswordAuthentication no" -o "StrictHostKeyChecking no"';
    ssh_options = [scp_options ' -x'];
    
    if isfield(Settings.SSH,'host_list')
        host = FindFreeSSH(Settings.SSH.host_list);
        if ~isempty(host)
            Settings.SSH.host = host;
        else
            error('openEMS:RunOpenEMS', 'unable to find host, abort openEMS');
        end
    end
    
    % create a tmp working dir
 	[status, result] = unix(['ssh ' ssh_options ' ' Settings.SSH.host ' "mktemp -d /tmp/openEMS_XXXXXXXXXXXX"']);
    if (status~=0)
        disp(result);
        error('openEMS:RunOpenEMS','mktemp failed to create tmp directory!');
    end
    ssh_work_path = strtrim(result); %remove tailing \n
    
    disp(['Running remote openEMS on ' Settings.SSH.host ' at working dir: ' ssh_work_path]);
    
    %copy openEMS all simulation files to the ssh host
    [stat, res] = unix(['scp ' scp_options ' * ' Settings.SSH.host ':' ssh_work_path '/']);
    if (stat~=0)
        disp(res);
        error('openEMS:RunOpenEMS','scp failed!');
    end

    %run openEMS (with log file if requested)
    if isfield(Settings,'LogFile')
        append_unix = [' 2>&1 | tee ' Settings.LogFile];
    else
        append_unix = [];
    end
	status = unix(['ssh ' ssh_options ' ' Settings.SSH.host ' "cd ' ssh_work_path ' && ' Settings.SSH.bin ' ' Sim_File ' ' opts '"' append_unix]);
    if (status~=0)
        disp(result);
        error('openEMS:RunOpenEMS','ssh openEMS failed!');
    end

    disp( 'Remote simulation done... copying back results and cleaning up...' );

    %copy back all results
    [stat, res] = unix(['scp -r ' scp_options ' ' Settings.SSH.host ':' ssh_work_path '/* ''' pwd '''/']);
    if (stat~=0);
        disp(res);
        error('openEMS:RunOpenEMS','scp failed!');
    end 
    
    %cleanup
    [stat, res] = unix(['ssh ' ssh_options ' ' Settings.SSH.host ' rm -r ' ssh_work_path]);
    if (stat~=0);
        disp(res);
        warning('openEMS:RunOpenEMS','remote cleanup failed!');
    end       
else
    args = [Sim_File ' ' opts];
    if isfield(Settings,'LogFile') && isfield(Settings,'Silent')
        invoke_openEMS(args,Settings.LogFile,Settings.Silent);
    elseif isfield(Settings,'LogFile')
        invoke_openEMS(args,Settings.LogFile);
    elseif isfield(Settings,'Silent')
        invoke_openEMS(args,[],Settings.Silent);
    else
        invoke_openEMS(args);
    end
end
 
cd(savePath);
return
