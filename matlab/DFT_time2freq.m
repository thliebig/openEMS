function f_val = DFT_time2freq( t, val, freq )
% f_val = DFT_time2freq( t, val, freq )
%
% computes the DFT at the given frequencies
% f_val: single-sided spectrum
%
% example:
%   t=linspace(0,1,100);
%   t_val=0.9*sin(2*pi*3*t); % sine wave; amplitude 0.9; frequency 3 Hz
%   f=linspace(1,5,101);
%   f_val=DFT_time2freq( t, t_val, f );
%   interp1(f,abs(f_val),3)
% ans = 0.8910
%   plot( t, t_val )
%   plot( f, abs(f_val) )
%
% See also FFT_time2freq

if numel(t) ~= numel(val)
    error 'numel(t) ~= numel(val)'
end

dt = t(2)-t(1);

f_val = zeros(1,numel(freq));
for f_idx=1:numel(freq)
    f_val(f_idx) = sum( val .* exp( -1i * 2*pi*freq(f_idx) * t ) );
end

% scaling for pulse-like signals
% as in the Fourier transform (https://en.wikipedia.org/wiki/Fourier_transform#Definition)
% the unit of the result is V/Hz, A/Hz, .../Hz
f_val = f_val * dt;

% scaling for periodic signals
% as in the Fourier series (https://en.wikipedia.org/wiki/Fourier_series#Definition)
% the unit of the result is V, A, ...
%f_val = f_val / length(t);
 
f_val = f_val * 2; % single-sided spectrum
