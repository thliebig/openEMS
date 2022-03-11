function host = FindFreeSSH(host_list, Settings, wait_time, command)
% function host = FindFreeSSH(host_list, Settings, wait_time, command)
% 
% Find a free ssh host not running openEMS
% 
% internal function used by RunOpenEMS
% 
% host_list: give a list of possible host
% 
% wait_time: wait x seconds after not finding a free host and rechecking
%            default: 600 seconds
% 
% command: unix command to check for free host (empty result --> free)
%          default: 'ps -e | grep openEMS'
%
% See also RunOpenEMS
%
% openEMS matlab interface
% -----------------------
% author: Thorsten Liebig

if (nargin<4)
    % command which should return an empty string if host is available
    command = 'ps -e | grep openEMS';
end

% 10 seconds ssh timeout
time_out = 10;

if (nargin<3)
    wait_time = 600;
end

if ~isunix
    ssh_command = [Settings.SSH.Putty.Path '/plink '];
    ssh_options = [' -i ' Settings.SSH.Putty.Key];
    command = ['"' command '"'];
else
    ssh_command = 'ssh';
    ssh_options = ['-o ConnectTimeout=' num2str(time_out)];
    command = ['''' command ''''];
end

if ischar(host_list)
    fid=fopen(host_list);
    if (fid==-1)
        error('FindFreeSSH: cannot open host file');
    end
    clear host_list;
    host_list = {};
    while 1
        line = fgetl(fid);
        if ischar(line)
            host_list{end+1} = line;
        else
            break;
        end
    end
    fclose(fid);
elseif ~iscell(host_list)
    error('FindFreeSSH: unknown host list format');
end

while 1
    for n = 1:numel(host_list)
        host = host_list{n};
        [status, result] = unix([ssh_command ' ' ssh_options ' ' host ' ' command ]);
        if (isempty(result) && status==1)
            disp(['FindFreeSSH:: found a free host: ' host ]);
            return
        elseif (~isempty(result) && status==0)
            disp(['FindFreeSSH:: ' host ' is busy running openEMS ... ' ]);
        else
            disp(['FindFreeSSH:: shh connection to ' host ' failed ... ' ]);
        end
    end
    
    host = '';
    
    if (wait_time<=0)
        warning('openEMS:FindFreeSSH',' unable to find a free host ');
        return
    end
    
    disp([' no free host found waiting for ' num2str(wait_time) ' seconds ... '])
    pause(wait_time)
end