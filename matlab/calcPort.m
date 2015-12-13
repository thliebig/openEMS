function [port] = calcPort( port, SimDir, f, varargin)
% [port] = calcPort( port, SimDir, f, varargin)
%
% Calculate:
%   - voltages and currents
%   - the propagation constant and the characteristic impedance (if applicable)
%
% The port has to be created by e.g. AddMSLPort(), AddLumpedPort() or AddCurvePort
%
% input:
%   port:       return value of AddMSLPort()
%   SimDir:     directory, where the simulation files are
%   f:          frequency vector for DFT
% 
% variable input:
%   'RefImpedance':  - use a given reference impedance to calculate inc and
%                      ref voltages and currents
%                    - default is given port or calculated line impedance
%   'RefPlaneShift': for transmission lines only, See also calcTLPort for
%                    more details
%   'SwitchDirection': 0/1, switch assumed direction of propagation
%   'SignalType':    'pulse' (default) or 'periodic'
%
% output: 
%   % output signals/values in time domain (TD):
%   port.ut.tot     total voltage (time-domain)
%   port.ut.time    voltage time vector
%   port.it.tot     total current (time-domain)
%   port.it.time    current time vector
%
%   % output signals/values in frequency domain (FD):
%   port.f                  the given frequency fector
%   port.uf.tot/inc/ref     total, incoming and reflected voltage
%   port.if.tot/inc/ref     total, incoming and reflected current
%   port.ZL_ref             used refernce impedance
%
%   port.P_inc              incoming power
%   port.P_ref              reflected power
%   port.P_acc              accepted power (incoming minus reflected,
%                           may be negative for passive ports)
%
%   if port is a transmission line port:
%   port.beta:              propagation constant
%   port.ZL:                characteristic line impedance
%
% example:
%   port = calcPort(port, Sim_Path, f, 'RefImpedance', 50);
%
% openEMS matlab interface
% -----------------------
% (C) 2012 Thorsten Liebig <thorsten.liebig@gmx.de>
%
% See also AddMSLPort, AddLumpedPort, AddCurvePort, calcTLPort, calcLumpedPort

if (iscell(port))
    for n=1:numel(port)
        port{n}=calcPort(port{n}, SimDir, f, varargin{:});
    end
    return;
end

if isempty(port)
    return;
end

if (strcmpi(port.type,'MSL') || strcmpi(port.type,'Coaxial') || strcmpi(port.type,'StripLine') || strcmpi(port.type,'CPW'))
    port = calcTLPort( port, SimDir, f, varargin{:});
elseif strcmpi(port.type,'WaveGuide')
    port = calcWGPort( port, SimDir, f, varargin{:});
elseif (strcmpi(port.type,'Lumped') || strcmpi(port.type,'Curve'))
    port = calcLumpedPort( port, SimDir, f, varargin{:});
else
    error 'unknown port type'
end

% calc some more port parameter
% incoming power
port.P_inc = 0.5*real(port.uf.inc.*conj(port.if.inc));
% reflected power
port.P_ref = 0.5*real(port.uf.ref.*conj(port.if.ref));
% accepted power (incoming - reflected)
port.P_acc = 0.5*real(port.uf.tot.*conj(port.if.tot));
