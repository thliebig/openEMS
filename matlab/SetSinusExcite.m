function FDTD = SetSinusExcite(FDTD,f0)
% function FDTD = SetSinusExcite(FDTD,f0)
%
% see also SetGaussExcite SetCustomExcite
%
% e.g FDTD = SetSinusExcite(FDTD,1e9);
%
% openEMS matlab interface
% -----------------------
% author: Thorsten Liebig

FDTD.Excitation.ATTRIBUTE.Type=1;
FDTD.Excitation.ATTRIBUTE.f0=f0;
FDTD.ATTRIBUTE.f_max=f0;
