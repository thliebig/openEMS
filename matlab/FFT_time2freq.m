function [f,val] = FFT_time2freq( t, val )
% [f,val] = FFT_time2freq( t, val )
%
% Note: This function can only be used for pulse signals
%
% See also DFT_time2freq

dt=t(2)-t(1); % timestep
L=numel(val); % signal length
NFFT = 2^nextpow2(L); % next power of 2 (makes fft fast)
%very fine freq resolution... NFFT = NFFT+100000;
val = fft( val, NFFT)*dt;
f = 1/(2*dt) * linspace(0,1,NFFT/2+1);

val = 2*val(1:NFFT/2+1); % single-sided spectrum

%correct phase for time-shifted signals
val = val .* exp(-1j*2*pi*f * t(1));
