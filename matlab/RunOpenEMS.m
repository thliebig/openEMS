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
% optinal:  
% (ssh only on unix with working ssh client)
% Settings.SSH.host = '<hostname or ip>'
% Settings.SSH.bin = '<path_to_openEMS>/openEMS.sh'
%
% Settings.LogFile = 'openEMS.log'
%
% RunOpenEMS(Sim_Path,Sim_File,opts,Settings)
%
% See also WriteOpenEMS
%
% openEMS matlab interface
% -----------------------
% author: Thorsten Liebig

if nargin < 3
    error 'specify the Sim_Path and Sim_file to simulate'
end

if (nargin<4)
    Settings = [];
end


savePath = pwd;
cd(Sim_Path);
    
if (isfield(Settings,'SSH') && isunix)
    ssh_work_path = ['openEMS_' int2str(randi([1e8 9e8],1))];

    disp(['Running remote openEMS on ' Settings.SSH.host ' at working dir: ' ssh_work_path]);

 	[status, result] = unix(['ssh ' Settings.SSH.host ' "mkdir /tmp/' ssh_work_path '"']);
    if (status~=0)
        disp(result);
        error('openEMS:RunOpenEMS','mkdir failed!');
    end
    
    [stat, res] = unix(['scp ' Sim_File ' ' Settings.SSH.host ':/tmp/' ssh_work_path '/' Sim_File]);
    if (stat~=0)
        disp(res);
        error('openEMS:RunOpenEMS','scp failed!');
    end

    if isfield(Settings,'LogFile')
        append_unix = [' 2>&1 | tee ' Settings.LogFile];
    else
        append_unix = [];
    end
	status = unix(['ssh ' Settings.SSH.host ' "cd /tmp/' ssh_work_path ' && ' Settings.SSH.bin ' ' Sim_File ' ' opts '"' append_unix])
    if (status~=0)
        disp(result);
        error('openEMS:RunOpenEMS','ssh openEMS failed!');
    end

    disp(['Remote simulation done... copying back results and cleaning up...']);

    [stat, res] = unix(['scp -r ' Settings.SSH.host ':/tmp/' ssh_work_path '/* ' pwd '/']);
    if (stat~=0);
        disp(res);
        error('openEMS:RunOpenEMS','scp failed!');
    end 
    
    [stat, res] = unix(['ssh ' Settings.SSH.host ' rm -r /tmp/' ssh_work_path]);
    if (stat~=0);
        disp(res);
        error('openEMS:RunOpenEMS','remote cleanup failed!');
    end       
else
    args = [Sim_File ' ' opts];
    if isfield(Settings,'LogFile')
        invoke_openEMS(args,Settings.LogFile);
    else
        invoke_openEMS(args);
    end
end
 
cd(savePath);
return
