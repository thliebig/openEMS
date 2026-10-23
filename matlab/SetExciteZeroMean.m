function FDTD = SetExciteZeroMean(FDTD,val)
% function FDTD = SetExciteZeroMean(FDTD,val);
%
% Force the excitation waveform to have zero time-integral.
%
% A soft E excitation deposits net charge proportional to the integral of
% its waveform. Inside a structure closed by PEC that leftover is a static
% field nothing can remove, and it stops the energy criterion from ever
% tripping. Set this when the excitation band reaches down toward DC.
%
% The correction subtracts a Hann-shaped bump scaled to cancel the integral
% exactly, applied to both the voltage and current signals. Off by default.
% It may be called before or after the Set*Excite function.
%
% - val : 1 (or true, default) to enable the correction, 0 to disable it
%
% see also SetGaussExcite SetCustomExcite
%
% e.g
%
%     FDTD = SetGaussExcite(FDTD,1.55e9,1.45e9);
%     FDTD = SetExciteZeroMean(FDTD,1);
%
% openEMS matlab interface
% -----------------------

if (nargin<2)
    val = 1;
end
FDTD.Excitation.ATTRIBUTE.ZeroMean=double(val~=0);
