function [port] = calcWGPort( port, SimDir, f, varargin)
% [port] = calcWGPort( port, SimDir, f, varargin)
%
% Calculate voltages and currents, the propagation constant beta
% and the characteristic impedance ZL of the given waveguide port.
%
% The port has to be created by e.g. AddWaveGuidePort().
%
% input:
%   port:       return value of e.g. AddWaveGuidePort()
%   SimDir:     directory, where the simulation files are
%   f:          frequency vector for DFT
%
% variable input:
%   'RefImpedance':  - use a given reference impedance to calculate inc and
%                      ref voltages and currents
%                    - default is given port or calculated line impedance
%   'RefPlaneShift': - use a given reference plane shift from port beginning
%                      for a desired phase correction
%                    - default is the measurement plane at the end of the
%                      port
%                    - the plane shift has to be given in drawing units!
%   'RefractiveIndex': set a material refractive index
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
%   port.f                  the given frequency fector
%   port.uf.tot/inc/ref     total, incoming and reflected voltage
%   port.if.tot/inc/ref     total, incoming and reflected current
%   port.beta:              propagation constant
%   port.ZL:                characteristic line impedance
%   port.ZL_ref             used reference impedance
%
% example:
%   port{1} = calcWGPort( port{1}, Sim_Path, f, 'RefImpedance', 50);
%
% openEMS matlab interface
% -----------------------
% (C) 2013 Thorsten Liebig (thorsten.liebig@gmx.de)
%
% See also AddWaveGuidePort, calcPort

if (iscell(port))
    for n=1:numel(port)
        port{n}=calcWGPort(port{n}, SimDir, f, varargin{:});
    end
    return;
end

if (strcmpi(port.type,'WaveGuide')~=1)
    error('openEMS:calcWGPort','error, type is not a waveguide port');
end

%set defaults
ref_ZL = -1;
ref_shift = nan;
ref_index = 1;
switch_dir = 1;

UI_args = {};

for n=1:2:numel(varargin)
    if (strcmp(varargin{n},'RefPlaneShift')==1);
        ref_shift = varargin{n+1};
    elseif (strcmp(varargin{n},'RefImpedance')==1);
        ref_ZL = varargin{n+1};
    elseif (strcmp(varargin{n},'RefractiveIndex')==1);
        ref_index = varargin{n+1};
    elseif (strcmpi(varargin{n},'SwitchDirection')==1);
        if (varargin{n+1})
            switch_dir = -1;
        end
    else
        UI_args(end+1) = varargin(n);
        UI_args(end+1) = varargin(n+1);
    end
end

% read time domain data
U = ReadUI( port.U_filename, SimDir, f, UI_args{:} );
I = ReadUI( port.I_filename, SimDir, f, UI_args{:} );

% store the original frequency domain waveforms
u_f = U.FD{1}.val;
i_f = I.FD{1}.val * switch_dir;

% time domain signal
port.ut.time  = U.TD{1}.t;
port.ut.tot = U.TD{1}.val;

port.it.time  = I.TD{1}.t;
port.it.tot = switch_dir*I.TD{1}.val;


physical_constants
k = 2*pi*f/C0*ref_index;
fc = C0*port.kc/2/pi/ref_index;
port.beta = sqrt(k.^2 - port.kc^2);
port.ZL = k * Z0 ./ port.beta;    %analytic waveguide impedance

% reference plane shift (lossless)
if ~isnan(ref_shift)
    % shift relative to the beginning of the waveguide
    ref_shift = ref_shift - port.measplanepos;
    ref_shift = ref_shift * port.drawingunit;

    % store the shifted frequency domain waveforms
    phase = real(beta)*ref_shift;
    u_f_shift = u_f .* cos(-phase) + 1i * i_f.*port.ZL .* sin(-phase);
    i_f_shift = i_f .* cos(-phase) + 1i * u_f./port.ZL .* sin(-phase);

    u_f = u_f_shift;
    i_f = i_f_shift;
end

if (ref_ZL < 0)
    ref_ZL = port.ZL;
end

port.ZL_ref = ref_ZL;

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
