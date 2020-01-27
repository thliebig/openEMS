function RunOpenEMS(Sim_Path, Sim_File, opts, Settings)
% function RunOpenEMS(Sim_Path, Sim_File, <opts, Settings>)
%
% Run an openEMS simulation.
%
% arguments:
% Sim_Path: specify the simulation folder (folder must exist!)
% Sim_File: xml-filename to simulate, created by WriteOpenEMS
%
% optional arguments
% 
% opts:  list of openEMS options
%        possible options:
%         --disable-dumps      Disable all field dumps for faster simulation
%         --debug-material     Dump material distribution to a vtk file for debugging
%         --debug-PEC          Dump metal distribution to a vtk file for debugging
%         --debug-operator     Dump operator to vtk file for debugging
%         --debug-boxes        Dump e.g. probe boxes to vtk file for debugging
%         --debug-CSX          Write CSX geometry file to debugCSX.xml
%         --engine=<type>      Choose engine type
%             --engine=fastest         fastest available engine (default)
%             --engine=basic           basic FDTD engine
%             --engine=sse             engine using sse vector extensions
%             --engine=sse-compressed  engine using compressed operator + sse vector extensions
%             --engine=MPI             engine using compressed operator + sse vector extensions + MPI parallel processing
%             --engine=multithreaded   engine using compressed operator + sse vector extensions + MPI + multithreading
%         --numThreads=<n>     Force use n threads for multithreaded engine
%         --no-simulation      only run preprocessing; do not simulate
%         --dump-statistics    dump simulation statistics to 'openEMS_run_stats.txt' and 'openEMS_stats.txt'
%
%          Additional global arguments
%         --showProbeDiscretization    Show probe discretization information
%         --nativeFieldDumps           Dump all fields using the native field components
%         -v,-vv,-vvv                  Set debug level: 1 to 3
%
%
% settings: list of Matlab settings
%       possible settings:
%           Settings.LogFile = 'openEMS.log'
%           Settings.Silent  = 0
%
%       additional remote simulation settings
%       Note: ssh only on unix with working ssh client or windows with putty client
%           openEMS Linux server or Windows with cygwin necessary
%       Settings.SSH.host = '<hostname or ip>'
%       Settings.SSH.bin = '<path_to_openEMS>/openEMS.sh'
%       ssh optional:
%       Settings.SSH.host_list = {'list','of','hosts'}; %searches for a free host
%       %on Windows needed additionally
%       Settings.SSH.Putty.Path = '<path_to>\putty';
%       Settings.SSH.Putty.Key = '<path_to>\putty_private_key.ppk';
%
%       MPI settings:
%       Settings.MPI.xxx --> help RunOpenEMS_MPI
% 
%
% example:
% %create CSX and FDTD
% WriteOpenEMS('/tmp/path_to_run_in/myfile.xml', FDTD, CSX)
% RunOpenEMS('/tmp/path_to_run_in','myfile.xml','-v')
%
% See also WriteOpenEMS FindFreeSSH InitCSX InitFDTD RunOpenEMS_MPI
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

ssh_command = 'ssh';
scp_command = 'scp';
scp_options = '';
ssh_options = '';

enable_ssh = 0;
enable_ssh = isfield(Settings,'SSH') && isunix;

if ~isunix
    enable_ssh = isfield(Settings,'SSH') && isfield(Settings.SSH,'Putty');
    if (enable_ssh)
        ssh_command = [Settings.SSH.Putty.Path '/plink '];
        ssh_options = [ssh_options ' -i ' Settings.SSH.Putty.Key];

        scp_command = [Settings.SSH.Putty.Path '/pscp '];
        scp_options = [scp_options ' -i ' Settings.SSH.Putty.Key];
    end
end

savePath = pwd;
cd(Sim_Path);

if (enable_ssh)
    scp_options = [scp_options ' -C'];
    ssh_options = [ssh_options ' -x -C'];

    % ssh options: no X forwarding; no password prompt (use pub keys!); no host checking
    if (isunix)
        ssh_options = [ssh_options ' -o "PasswordAuthentication no" -o "StrictHostKeyChecking no"'];
        scp_options = [scp_options ' -o "PasswordAuthentication no" -o "StrictHostKeyChecking no"'];
    end

    if isfield(Settings.SSH,'host_list')
        host = FindFreeSSH(Settings.SSH.host_list, Settings);
        if ~isempty(host)
            Settings.SSH.host = host;
        else
            error('openEMS:RunOpenEMS', 'unable to find host, abort openEMS');
        end
    end

    % create a tmp working dir
    [status, result] = system([ssh_command ' ' ssh_options ' ' Settings.SSH.host ' "mktemp -d /tmp/openEMS_XXXXXXXXXXXX"']);
    if (status~=0)
        disp(result);
        error('openEMS:RunOpenEMS','mktemp failed to create tmp directory!');
    end
    ssh_work_path = strtrim(result); %remove tailing \n

    disp(['Running remote openEMS on ' Settings.SSH.host ' at working dir: ' ssh_work_path]);

    %copy openEMS all simulation files to the ssh host
    [stat, res] = system([scp_command ' ' scp_options ' * ' Settings.SSH.host ':' ssh_work_path '/']);
    if (stat~=0)
        disp(res);
        error('openEMS:RunOpenEMS','scp failed!');
    end

    %run openEMS (with log file if requested)
    if isfield(Settings,'LogFile') && isunix
        append_unix = [' 2>&1 | tee ' Settings.LogFile];
    else
        append_unix = [];
    end
    status = system([ssh_command ' ' ssh_options ' ' Settings.SSH.host ' "cd ' ssh_work_path ' && ' Settings.SSH.bin ' ' Sim_File ' ' opts '"' append_unix]);
    if (status~=0)
        disp(result);
        error('openEMS:RunOpenEMS','ssh openEMS failed!');
    end

    disp( 'Remote simulation done... copying back results and cleaning up...' );

    %copy back all results
    [stat, res] = system([scp_command ' -r ' scp_options ' ' Settings.SSH.host ':' ssh_work_path '/* ./']);
    if (stat~=0);
        disp(res);
        error('openEMS:RunOpenEMS','scp failed!');
    end

    %cleanup
    [stat, res] = system([ssh_command ' ' ssh_options ' ' Settings.SSH.host ' rm -r ' ssh_work_path]);
    if (stat~=0);
        disp(res);
        warning('openEMS:RunOpenEMS','remote cleanup failed!');
    end
else
    args = ['"' Sim_File '" ' opts]
    if isfield(Settings,'LogFile') && isfield(Settings,'Silent')
        invoke_openEMS(args,['"' Settings.LogFile '"'],Settings.Silent);
    elseif isfield(Settings,'LogFile')
        invoke_openEMS(args,['"' Settings.LogFile '"']);
    elseif isfield(Settings,'Silent')
        invoke_openEMS(args,[],Settings.Silent);
    else
        invoke_openEMS(args);
    end
end

cd(savePath);
return
