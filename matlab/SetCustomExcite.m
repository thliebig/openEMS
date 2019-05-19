function FDTD = SetCustomExcite(FDTD,f0,funcStr)
% function FDTD = SetCustomExcite(FDTD,f0,funcStr)
%
% f0 : nyquist rate
% funcStr : string describing the excitation function e(t)
%
% see also SetSinusExcite SetGaussExcite
%
% e.g for a ramped sinus excite...
% T = 1/f0;
% FDTD = SetCustomExcite(FDTD,1e9,..
% [ '(1-exp(-1*(t/' num2str(T) ')^2) ) * sin(2*pi*' num2str(f0) '*t)' ]);
%
% openEMS matlab interface
% -----------------------
% author: Thorsten Liebig

FDTD.Excitation.ATTRIBUTE.Type=10;
FDTD.Excitation.ATTRIBUTE.f0=f0;
FDTD.Excitation.ATTRIBUTE.Function=funcStr;
FDTD.ATTRIBUTE.f_max=f0;
