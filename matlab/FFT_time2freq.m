function [f,val] = FFT_time2freq( t, val )

dt=t(2)-t(1);
val = [val zeros(1,5000)];
L=numel(val);
f = (0:L-1)/L/dt;
f = f(1:floor(L/2));
val = 2*fft(val)/L;
val = val(1:floor(L/2));
