function [f,val] = FFT_time2freq( t, val )

dt=t(2)-t(1); % timestep
L=numel(val); % signal length
NFFT = 2^nextpow2(L); % next power of 2 (makes fft fast)
%very fine freq resolution... NFFT = NFFT+100000;
val = fft( val, NFFT)/L;
f = 1/(2*dt) * linspace(0,1,NFFT/2+1);

val = 2*val; % single-sided spectrum
