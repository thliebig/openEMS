function [S11,beta,ZL] = calcPort( portstruct, SimDir, f, ref_shift )
%[S11,beta,ZL] = calcMSLPort( portstruct, SimDir, [f], [ref_shift] )
%
% Calculate the reflection coefficient S11, the propagation constant beta
% of the MSL-port and the characteristic impedance ZL of the MSL-port.
% The port is to be created by AddMSLPort().
%
% input:
%   portstruct: return value of AddMSLPort()
%   SimDir:     directory, where the simulation files are
%   f:          (optional) frequency vector for DFT
%   ref_shift:  (optional) reference plane shift measured from start of port (in drawing units)
%
% output:
%   S11:  reflection coefficient (normalized to ZL)
%   beta: propagation constant
%   ZL:   characteristic line impedance
%
% reference: W. K. Gwarek, "A Differential Method of Reflection Coefficient Extraction From FDTD Simulations",
%            IEEE Microwave and Guided Wave Letters, Vol. 6, No. 5, May 1996
%
% openEMS matlab interface
% -----------------------
% (C) 2010 Sebastian Held <sebastian.held@uni-due.de>
% See also AddMSLPort

%DEBUG
% save('/tmp/test.mat', 'portstruct', 'SimDir', 'f', 'nargin' )
% load('/tmp/test.mat')

% check
if portstruct.v_delta(1) ~= portstruct.v_delta(2)
	warning( 'openEMS:calcPort:mesh', 'mesh is not equidistant; expect degraded accuracy' );
end

if nargin < 3
    f = [];
end

% read time domain data
filename = ['port_ut' num2str(portstruct.nr)];
U = ReadUI( {[filename 'A'],[filename 'B'],[filename 'C']}, SimDir, f );
filename = ['port_it' num2str(portstruct.nr)];
I = ReadUI( {[filename 'A'],[filename 'B']}, SimDir, f );

f = U.FD{2}.f;
Et = U.FD{2}.val;
dEt = (U.FD{3}.val - U.FD{1}.val) / (sum(abs(portstruct.v_delta(1:2))) * portstruct.drawingunit);
Ht = (I.FD{1}.val + I.FD{2}.val)/2; % space averaging: Ht is now defined at the same pos as Et
dHt = (I.FD{2}.val - I.FD{1}.val) / (abs(portstruct.i_delta(1)) * portstruct.drawingunit);

beta = sqrt( - dEt .* dHt ./ (Ht .* Et) );
beta(real(beta) < 0) = -beta(real(beta) < 0); % determine correct sign (unlike the paper)

% determine S11
A = sqrt( Et .* dHt ./ (Ht .* dEt) );
A(imag(A) > 0) = -A(imag(A) > 0); % determine correct sign (unlike the paper)
S11 = (A - 1) ./ (A + 1);

% determine S11_corrected
delta_e = sum(portstruct.v_delta(1:2))/2 * portstruct.drawingunit;
delta_h = portstruct.i_delta(1) * portstruct.drawingunit;
S11_corrected = sqrt( Et .* (dHt ./ (sin(beta.*delta_h*.5)/(beta*delta_h*.5))) ./ ((Ht ./ cos(beta*delta_h*.5)) .* (dEt ./ (sin(beta*delta_e)./(beta*delta_e)))));
S11_corrected(imag(S11_corrected) > 0) = -S11_corrected(imag(S11_corrected) > 0); % determine correct sign (unlike the paper)
S11_corrected = (S11_corrected-1) ./ (S11_corrected+1);

% my own solution...
temp = sqrt(-dHt .* dEt ./ (Ht .* Et));
S11 = (-1i * dEt + Et .* temp) ./ (Et .* temp + 1i * dEt); % solution 1
% S11 = (-1i * dEt - Et .* temp) ./ (-Et .* temp + 1i * dEt); % solution 2

% % determine ZL
% Et_forward = Et ./ (1 + S11);
% Ht_forward = Ht ./ (1 - S11);
% ZL = Et_forward ./ Ht_forward;
% 
% % determine ZL_corrected
% Et_forward_corrected = Et ./ (1 + S11_corrected);
% Ht_forward_corrected = Ht ./ (1 - S11_corrected);
% ZL_corrected = Et_forward_corrected ./ Ht_forward_corrected;

% determine ZL
ZL = sqrt(Et .* dEt ./ (Ht .* dHt));

% reference plane shift (lossless)
if (nargin > 3)
    % renormalize the shift to the measurement plane
    ref_shift = ref_shift - portstruct.measplanepos;
    ref_shift = ref_shift * portstruct.drawingunit;
    S11 = S11 .* exp(2i*real(beta)*ref_shift);
    S11_corrected = S11_corrected .* exp(2i*real(beta)*ref_shift);
end
