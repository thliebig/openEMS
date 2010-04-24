function WriteOpenEMS(filename, FDTD, CSX)
% function WriteOpenEMS(filename, FDTD, CSX)
% 
% openEMS matlab interface
% -----------------------
% author: Thorsten Liebig

openEMS.FDTD = FDTD;
openEMS.ContinuousStructure = CSX;
struct_2_xml(filename,openEMS,'openEMS');