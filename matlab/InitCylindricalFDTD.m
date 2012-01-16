function FDTD = InitCylindricalFDTD(NrTS, endCrit, varargin)
% function FDTD = InitCylindricalFDTD(NrTS, endCrit, varargin)
%
% see also InitFDTD
%
% e.g FDTD = InitCylindricalFDTD(5e5,1e-6,'OverSampling',10)
% 
% WARNING: This function is depreciated, use InitFDTD with 'CoordSystem',1
%          e.g.: InitFDTD(5e5,1e-6,'OverSampling',10, 'CoordSystem',1)
%
% openEMS matlab interface
% -----------------------
% author: Thorsten Liebig

warning('InitCylindricalFDTD: This function is depreciated, use InitFDTD with ''CoordSystem'',1');

FDTD = InitFDTD(NrTS, endCrit, varargin{:});

FDTD.ATTRIBUTE.CylinderCoords=1;


