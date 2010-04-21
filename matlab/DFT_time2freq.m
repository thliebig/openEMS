function f_val = DFT_time2freq( t, val, freq )
% val = FFT_time2freq( t, val, freq )
%
% computes the DFT at the given frequencies

if numel(t) ~= numel(val)
    error 'numel(t) ~= numel(val)'
end

f_val = zeros(1,numel(freq));
for f_idx=1:numel(freq)
    f_val(f_idx) = sum( val .* exp( -1i * 2*pi*freq(f_idx) * t ) );
end
f_val = f_val / numel(t);

f_val = f_val * 2; % single-sided spectrum
