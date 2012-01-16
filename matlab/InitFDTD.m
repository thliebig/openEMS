function FDTD = InitFDTD(NrTS, endCrit, varargin)
% function FDTD = InitFDTD(NrTS, endCrit, varargin)
%
% possible arguments for useage with openEMS:
%   OverSampling:   nyquist oversampling of time domain dumps
%   CoordSystem:    choose coordinate system (0 Cartesian, 1 Cylindrical)
%
% e.g FDTD = InitFDTD(5e5,1e-6,'OverSampling',10)
%
% openEMS matlab interface
% -----------------------
% author: Thorsten Liebig

FDTD.ATTRIBUTE.NumberOfTimesteps=NrTS;
FDTD.ATTRIBUTE.endCriteria=endCrit;

for n=1:numel(varargin)/2
    if strcmp(varargin{2*n-1},'CoordSystem')==1
        FDTD.ATTRIBUTE.CylinderCoords=varargin{2*n}==1;
    else
        FDTD.ATTRIBUTE.(varargin{2*n-1})=varargin{2*n};
    end
end

