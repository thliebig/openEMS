function f_val = DFT_time2freq( t, val, freq, signal_type )
% f_val = DFT_time2freq( t, val, freq, signal_type )
%
% computes the DFT at the given frequencies
%
% parameter:
% t  : time vector
% val: data vector
% freq: DFT frequency vector
% signal_type: 'pulse' (default), 'periodic'
%
% return values:
% f_val:  single-sided spectrum
%
% example:
%   t=linspace(0,1,100);
%   t_val=0.9*sin(2*pi*3*t); % sine wave; amplitude 0.9; frequency 3 Hz
%   f=linspace(1,5,101);
%   f_val=DFT_time2freq( t, t_val, f, 'periodic' );
%   interp1(f,abs(f_val),3)
% ans = 0.8910
%   plot( t, t_val )
%   plot( f, abs(f_val) )

if numel(t) ~= numel(val)
    error 'numel(t) ~= numel(val)'
end

if nargin<4
  signal_type = 'pulse';
end

f_val = zeros(1,numel(freq));
for f_idx=1:numel(freq)
    f_val(f_idx) = sum( val .* exp( -1i * 2*pi*freq(f_idx) * t ) );
end

if strcmpi(signal_type, 'pulse')
  dt = t(2)-t(1);
  f_val = f_val * dt;
elseif strcmpi(signal_type, 'periodic')
  f_val = f_val / length(t);
else
  error 'unknown signal type'
end
 
f_val = f_val * 2; % single-sided spectrum
