%close all;
clear all;
clc

tmp = load('../tmp/u1');

t = tmp(:,1);
u = tmp(:,2);

L=numel(t);

subplot(2,1,1);
plot(t,u);

dt=t(2)-t(1);

f = (1:L)/L/dt;
fu = fft(u)/L;
subplot(2,1,2);
plot(f(1:L/2),abs(fu(1:L/2)));

