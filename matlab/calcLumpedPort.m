function [port] = calcLumpedPort( port, SimDir, f, varargin)
% [port] = calcLumpedPort( port, SimDir, f, varargin)
%
% Calculate voltages and currents of given lumped port.
%
% The port has to be created by e.g. AddLumpedPort().
%
% input:
%   port:       return value of e.g. AddMSLPort()
%   SimDir:     directory, where the simulation files are
%   f:          frequency vector for DFT
%
% variable input:
%   'RefImpedance':  - use a given reference impedance to calculate inc and
%                      ref voltages and currents
%                    - default is given port or calculated line impedance
%   'SwitchDirection': 0/1, switch assumed direction of propagation
%
% output:
%   % output signals/values in time domain (TD):
%   port.ut.tot     total voltage (time-domain)
%   port.ut.time    voltage time vector
%   port.it.tot     total current (time-domain)
%   port.it.time    current time vector
%
%   % output signals/values in frequency domain (FD):
%   port.f                  the given frequency factor
%   port.uf.tot/inc/ref     total, incoming and reflected voltage
%   port.if.tot/inc/ref     total, incoming and reflected current
%
% example:
%   port{1} = calcLumpedPort( port{1}, Sim_Path, f, 'RefImpedance', 50);
%
% openEMS matlab interface
% -----------------------
% (C) 2012 Thorsten Liebig <thorsten.liebig@gmx.de>
%
% See also AddLumpedPort, calcPort

if (iscell(port))
    for n=1:numel(port)
        port{n}=calcLumpedPort(port{n}, SimDir, f, varargin{:});
    end
    return;
end

if (strcmpi(port.type,'Lumped')~=1 && strcmpi(port.type,'Curve')~=1)
    error('openEMS:calcLumpedPort','error, type is not a lumped port');
end


%% read optional arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%set defaults
ref_ZL = port.Feed_R;
switch_dir = 1;

UI_args = {};

for n=1:2:numel(varargin)
    if (strcmp(varargin{n},'RefImpedance')==1);
        ref_ZL = varargin{n+1};
    elseif (strcmpi(varargin{n},'SwitchDirection')==1);
        if (varargin{n+1})
            switch_dir = -1;
        end
    else
        UI_args(end+1) = varargin(n);
        UI_args(end+1) = varargin(n+1);
    end
end

port.ZL_ref = ref_ZL;

% read time domain data
U = ReadUI( port.U_filename, SimDir, f, UI_args{:} );
I = ReadUI( port.I_filename, SimDir, f, UI_args{:} );

% store the original frequency domain waveforms
u_f = U.FD{1}.val;
i_f = switch_dir*I.FD{1}.val;

port.ut.time  = U.TD{1}.t;
port.ut.tot = U.TD{1}.val;

port.it.time  = I.TD{1}.t;
port.it.tot = switch_dir*I.TD{1}.val;

port.Zin = u_f./i_f;

port.f = f;
uf_inc = 0.5 * ( u_f + i_f .* ref_ZL );
if_inc = 0.5 * ( i_f + u_f ./ ref_ZL );

uf_ref = u_f - uf_inc;
if_ref = if_inc - i_f; 

port.uf.tot = u_f;
port.uf.inc = uf_inc;
port.uf.ref = uf_ref;

port.if.tot = i_f;
port.if.inc = if_inc;
port.if.ref = if_ref;

port.raw.U = U;
port.raw.I = I;
