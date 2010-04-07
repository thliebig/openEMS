function FDTD = SetCustomExcite(FDTD,f0,funcStr);

FDTD.Excitation.ATTRIBUTE.Type=10;
FDTD.Excitation.ATTRIBUTE.f0=f0;
FDTD.Excitation.ATTRIBUTE.Function=funcStr;