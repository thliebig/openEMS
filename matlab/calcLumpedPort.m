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
%
% output:
%   port.f                  the given frequency fector
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
n_conv_arg = 3; % number of conventional arguments

%set defaults
ref_ZL = port.Feed_R;

if (nargin>n_conv_arg)
    for n=1:2:(nargin-n_conv_arg)
        if (strcmp(varargin{n},'RefImpedance')==1);
            ref_ZL = varargin{n+1};
        end
    end
end

% read time domain data
filename = ['port_ut' num2str(port.nr)];
U = ReadUI(filename, SimDir, f );
filename = ['port_it' num2str(port.nr)];
I = ReadUI(filename, SimDir, f );

% store the original frequency domain waveforms
u_f = U.FD{1}.val;
i_f = I.FD{1}.val; % shift to same position as v

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
