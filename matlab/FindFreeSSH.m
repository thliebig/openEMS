function host = FindFreeSSH(host_list, wait_time, command)
% function host = FindFreeSSH(host_list, wait_time, command)
% 
% Find a free ssh host not running openEMS
% 
% internal function used by RunOpenEMS
% 
% host_list: give a list of possible host
% 
% wait_time: wait x seconds after not finding a free host and rechecking
%            default: 60 seconds
% 
% command: unix command to check for free host (empty result --> free)
%          default: 'ps -ewwo user,args | grep openEMS'
%
% See also RunOpenEMS
%
% openEMS matlab interface
% -----------------------
% author: Thorsten Liebig

if (nargin<3)
    % command which should return an empty string if host is available
    command = 'ps -e user,args | grep openEMS';
end

if (nargin<2)
    wait_time = 60;
end

while 1
    
    for n = 1:numel(host_list)
        host = host_list{n};
        [status, result] = unix(['ssh ' host ' ' command]);
        
        if isempty(result)
            disp(['FindFreeSSH:: found a free host: ' host ]);
            return
        else
            disp(['FindFreeSSH:: ' host ' is busy running openEMS ... ' ]);
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