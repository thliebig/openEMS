function WriteOpenEMS(filename, FDTD, CSX)
% function WriteOpenEMS(filename, FDTD, CSX)
%
% Write the FDTD and CSX structures to a file.
%
% example:
% CSX = InitCSX();
% FDTD = InitFDTD();
% WriteOpenEMS('test.xml',FDTD,CSX)
%
% See also InitFDTD InitCSX CSXGeomPlot
%
% openEMS matlab interface
% -----------------------
% author: Thorsten Liebig

openEMS.FDTD = FDTD;
openEMS.ContinuousStructure = CSX;
struct_2_xml(filename,openEMS,'openEMS');