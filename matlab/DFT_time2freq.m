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

f_val = f_val * dt;
 
f_val = f_val * 2; % single-sided spectrum
