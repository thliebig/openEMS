function WriteOpenEMS(filename, FDTD, CSX)

openEMS.FDTD = FDTD;
openEMS.ContinuousStructure = CSX;
struct_2_xml(filename,openEMS,'openEMS');