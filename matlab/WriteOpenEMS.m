function WriteOpenEMS(filename, FDTD, CSX)

openEMS.FDTD=FDTD;
openEMS.ContinuousStructure=CSX;

xml_write(filename,openEMS);