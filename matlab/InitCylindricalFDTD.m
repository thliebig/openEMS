function FDTD = InitCylindricalFDTD(NrTS, endCrit, varargin)
% function FDTD = InitCylindricalFDTD(NrTS, endCrit, varargin)
%
% see also InitFDTD
%
% e.g FDTD = InitCylindricalFDTD(5e5,1e-6,'OverSampling',10)
%
% openEMS matlab interface
% -----------------------
% author: Thorsten Liebig

FDTD = InitFDTD(NrTS, endCrit, varargin);

FDTD.ATTRIBUTE.CylinderCoords=1;


