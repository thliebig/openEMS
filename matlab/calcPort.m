function [port] = calcPort( port, SimDir, f, varargin)
% [port] = calcPort( port, SimDir, f, varargin)
%
% Calculate voltages and currents, the propagation constant beta
% and the characteristic impedance ZL of the given port.
% The port has to be created by e.g. AddMSLPort().
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
%   'RefPlaneShift': - use a given reference plane shift from port beginning
%                      for a desired phase correction
%                    - default is the measurement plane
%                    - the plane shift has to be given in drawing units!
%
% output: 
%   port.f                  the given frequency fector
%   port.uf.tot/inc/ref     total, incoming and reflected voltage
%   port.if.tot/inc/ref     total, incoming and reflected current
%   port.beta:              propagation constant
%   port.ZL:                characteristic line impedance
%
% example:
%   port{1} = calcPort( port{1}, Sim_Path, f, 'RefImpedance', 50);
%  or 
% 
% reference: W. K. Gwarek, "A Differential Method of Reflection Coefficient Extraction From FDTD Simulations",
%            IEEE Microwave and Guided Wave Letters, Vol. 6, No. 5, May 1996
%
% openEMS matlab interface
% -----------------------
% (C) 2010 Sebastian Held <sebastian.held@uni-due.de>
% See also AddMSLPort

%DEBUG
% save('/tmp/test.mat', 'port', 'SimDir', 'f', 'nargin' )
% load('/tmp/test.mat')

% check
if abs((port.v_delta(1) - port.v_delta(2)) / port.v_delta(1))>1e-6
	warning( 'openEMS:calcPort:mesh', 'mesh is not equidistant; expect degraded accuracy' );
end


%% read optional arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_conv_arg = 3; % number of conventional arguments

%set defaults
ref_ZL = 0;
ref_shift = nan;

if (nargin>n_conv_arg)
    for n=1:2:(nargin-n_conv_arg)
        if (strcmp(varargin{n},'RefPlaneShift')==1);
            ref_shift = varargin{n+1};
        end
        
        if (strcmp(varargin{n},'RefImpedance')==1);
            ref_ZL = varargin{n+1};
        end
    end
end

% read time domain data
filename = ['port_ut' num2str(port.nr)];
U = ReadUI( {[filename 'A'],[filename 'B'],[filename 'C']}, SimDir, f );
filename = ['port_it' num2str(port.nr)];
I = ReadUI( {[filename 'A'],[filename 'B']}, SimDir, f );

% store the original frequency domain waveforms
u_f = U.FD{2}.val;
i_f = (I.FD{1}.val + I.FD{2}.val) / 2; % shift to same position as v

f = U.FD{2}.f;
Et = U.FD{2}.val;
dEt = (U.FD{3}.val - U.FD{1}.val) / (sum(abs(port.v_delta(1:2))) * port.drawingunit);
Ht = (I.FD{1}.val + I.FD{2}.val)/2; % space averaging: Ht is now defined at the same pos as Et
dHt = (I.FD{2}.val - I.FD{1}.val) / (abs(port.i_delta(1)) * port.drawingunit);

beta = sqrt( - dEt .* dHt ./ (Ht .* Et) );
beta(real(beta) < 0) = -beta(real(beta) < 0); % determine correct sign (unlike the paper)

% determine ZL
ZL = sqrt(Et .* dEt ./ (Ht .* dHt));

% reference plane shift (lossless)
if ~isnan(ref_shift)
    ref_shift = ref_shift * port.LengthScale;
    % shift to the beginning of MSL
    ref_shift = ref_shift - port.measplanepos;
    ref_shift = ref_shift * port.drawingunit;
    
    % store the shifted frequency domain waveforms
    phase = real(beta)*ref_shift;
    U.FD{1}.val = u_f .* cos(-phase) + 1i * i_f.*ZL .* sin(-phase);
    I.FD{1}.val = i_f .* cos(-phase) + 1i * u_f./ZL .* sin(-phase);
    
    u_f = U.FD{1}.val;
    i_f = I.FD{1}.val;
end

if (ref_ZL == 0)
    if isfield(port,'Feed_R')
        ref_ZL = port.Feed_R;
    else
        ref_ZL = ZL;
    end
end

port.ZL =  ZL;
port.beta = beta;

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
