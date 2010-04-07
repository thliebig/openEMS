function FDTD = InitFDTD(NrTS, endCrit, varargin)
% function FDTD = InitFDTD(NrTS, endCrit, varargin)
%
% possible arguments for useage with openEMS:
%   OverSampling:   nyquist oversampling of time domain dumps
%
% e.g FDTD = InitFDTD(5e5,1e-6,'OverSampling',10)
%
% openEMS matlab interface
% -----------------------
% author: Thorsten Liebig

FDTD.ATTRIBUTE.NumberOfTimesteps=NrTS;
FDTD.ATTRIBUTE.endCriteria=endCrit;

for n=1:(nargin-2)/2
    FDTD.ATTRIBUTE.(varargin{2*n-1})=varargin{2*n};
end

