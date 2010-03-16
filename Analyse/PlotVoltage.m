%close all;
clear all;
clc

fmax = 50e6;

figure(1);
tmpu = load('../tmp/u1');
tmpi = load('../tmp/i1');

t = tmpu(:,1);
u = tmpu(:,2);

subplot(2,2,1);
title('u_1 TD');
plot(t,u);
xlabel('t \rightarrow');
ylabel('ut_1 \rightarrow');
grid on;

dt=t(2)-t(1);
u= [u ; zeros(5000,1)];
L=numel(u);
t = (1:L)*dt;

f = (0:L-1)/L/dt;
fu = fft(u)/L;
subplot(2,2,2);
title('u_1 FD');
plot(f(1:L/2),abs(fu(1:L/2)));
xlabel('f \rightarrow');
ylabel('|uf_1| \rightarrow');
grid on;


t = tmpi(:,1);
i = tmpi(:,2);

subplot(2,2,3);
title('i_1 TD');
plot(t,i);
xlabel('t \rightarrow');
ylabel('it_1 \rightarrow');
grid on;

dt=t(2)-t(1);
i = [i; zeros(5000,1)];
L=numel(i);
t = (1:L)*dt;
f = (0:L-1)/L/dt;

fi = fft(i)/L;
subplot(2,2,4);
title('i_1 FD');
plot(f(1:L/2),abs(fi(1:L/2)));
xlabel('f \rightarrow');
ylabel('|if_1| \rightarrow');
grid on;

figure(2);
subplot(2,1,1);
plot(f,real(fu./fi));
xlim([0 fmax]);
grid on;
subplot(2,1,2);
plot(f,imag(fu./fi));
xlim([0 fmax]);
grid on;


