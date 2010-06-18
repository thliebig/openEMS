function Y = calc_ypar( f, ports, Sim_Path_Prefix )
% Y = calc_ypar( f, ports, Sim_Path_Prefix )
%
% f: frequency vector (Hz)
% ports: cell array of ports (see AddMSLPort() and AddLumpedPort())
% Sim_Path_Prefix: prefix of the simulation dirs (will be postfixed by
% excitation port number)
%
% This function calculates the Y-matrix representation of the ports
%
% It is assumed that each port (inside ports) is excited and the
% corresponding simulation was carried out at Sim_Path + portnr (e.g. for
% port 2: '/tmp/sim2')
%
% Sebastian Held <sebastian.held@uni-due.de>
% Jun 9 2010
%
% See also AddMSLPort AddLumpedPort

% sanitize input arguments
f = reshape(f,1,[]); % make it a row vector

% prepare result matrix
maxportnr = max( cellfun(@(x) x.nr, ports) );
Y = ones(maxportnr,maxportnr,numel(f)) * NaN;
U = ones(maxportnr,maxportnr,numel(f)) * NaN;
I = ones(maxportnr,maxportnr,numel(f)) * NaN;

% read time domain simulation results
for runnr = 1:numel(ports)
    Sim_Path = [Sim_Path_Prefix num2str(ports{runnr}.nr)];
    for pnr = 1:numel(ports)
        if isfield( ports{pnr}, 'v_delta' )
            % this is an MSLPort
            temp_U = ReadUI( ['port_ut' num2str(ports{pnr}.nr) 'B'], Sim_Path );
            temp = ReadUI( {['port_it' num2str(ports{pnr}.nr) 'A'],['port_it' num2str(ports{pnr}.nr) 'B']}, Sim_Path );
            temp_I.TD{1}.t = temp.TD{1}.t;
            temp_I.TD{1}.val = (temp.TD{1}.val + temp.TD{2}.val) / 2; % space averaging
        else
            % this is a lumped port
            temp_U = ReadUI( ['port_ut' num2str(ports{pnr}.nr)], Sim_Path );
            temp_I = ReadUI( ['port_it' num2str(ports{pnr}.nr)], Sim_Path );

%             % correct the orientation of the probes (FIXME to be done inside
%             % openEMS)
%             temp_U.TD{1}.val = temp_U.TD{1}.val * (-ports{pnr}.direction);
        end

%         % correct the orientation of the probes (FIXME to be done inside
%         % openEMS)
%         temp_I.TD{1}.val = temp_I.TD{1}.val * ports{pnr}.direction;
%         if runnr == 5 % DEBUG
%             temp_I.TD{1}.val = temp_I.TD{1}.val * -1;
%         end
        
        % time domain -> frequency domain
        U(ports{pnr}.nr,ports{runnr}.nr,:) = DFT_time2freq( temp_U.TD{1}.t, temp_U.TD{1}.val, f );
        I(ports{pnr}.nr,ports{runnr}.nr,:) = DFT_time2freq( temp_I.TD{1}.t, temp_I.TD{1}.val, f );

        % compensate H-field time advance
        delta_t_2 = temp_I.TD{1}.t(1) - temp_U.TD{1}.t(1); % half time-step (s) 
        I(ports{pnr}.nr,ports{runnr}.nr,:) = squeeze(I(ports{pnr}.nr,ports{runnr}.nr,:)).' .* exp(-1i*2*pi*f*delta_t_2);
    end
end

% calc Y-parameters
for a=1:numel(f)
    Y(:,:,a) = I(:,:,a) / U(:,:,a);
end
