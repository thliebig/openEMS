function FDTD = InitFDTD(NrTS, endCrit, varargin)
% function FDTD = InitFDTD(NrTS, endCrit, varargin)
%
% Inititalize the FDTD data-structure.
%
% optional arguments:
%   NrTS:       max. number of timesteps to simulate (e.g. default=1e9)
%   endCrit:    end criteria, e.g. 1e-5, simulations stops if energy has
%               decayed by this value (<1e-4 is recommended, default=1e-5)
%
% optional field arguments for usage with openEMS:
%   OverSampling:   nyquist oversampling of time domain dumps
%   CoordSystem:    choose coordinate system (0 Cartesian, 1 Cylindrical)
%
% examples:
%   %default init with 1e9 max. timesteps and -50dB end-criteria
%   FDTD = InitFDTD();
%
%   %init with 1e6 max. timesteps and -60dB end-criteria
%   FDTD = InitFDTD(1e6, 1e-6);
%
%   %cylindrical FDTD simulation
%   FDTD = InitFDTD(1e6, 1e-6, 'CoordSystem', 1);
%
%   See also InitCSX
%
% openEMS matlab interface
% -----------------------
% author: Thorsten Liebig (c) 2010-2012

if (nargin<1)
   FDTD.ATTRIBUTE.NumberOfTimesteps=1e9;
else
    FDTD.ATTRIBUTE.NumberOfTimesteps=NrTS;
end
if (nargin<2)
    FDTD.ATTRIBUTE.endCriteria = 1e-5;
else
    FDTD.ATTRIBUTE.endCriteria=endCrit;
end

for n=1:numel(varargin)/2
    if strcmp(varargin{2*n-1},'CoordSystem')==1
        FDTD.ATTRIBUTE.CylinderCoords=varargin{2*n}==1;
    else
        FDTD.ATTRIBUTE.(varargin{2*n-1})=varargin{2*n};
    end
end

