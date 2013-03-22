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
%
% output: 
%   port.f                  the given frequency fector
%   port.uf.tot/inc/ref     total, incoming and reflected voltage
%   port.if.tot/inc/ref     total, incoming and reflected current
%   port.ZL_ref             used refernce impedance
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

if strcmpi(port.type,'MSL')
    port = calcTLPort( port, SimDir, f, varargin{:});
    return
elseif (strcmpi(port.type,'Lumped') || strcmpi(port.type,'Curve'))
    port = calcLumpedPort( port, SimDir, f, varargin{:});
    return
else
    error 'unknown port type'
end

