function FDTD = SetBoundaryCond(FDTD,BC)

FDTD.BoundaryCond.ATTRIBUTE.xmin=BC(1);
FDTD.BoundaryCond.ATTRIBUTE.xmax=BC(2);
FDTD.BoundaryCond.ATTRIBUTE.ymin=BC(3);
FDTD.BoundaryCond.ATTRIBUTE.ymax=BC(4);
FDTD.BoundaryCond.ATTRIBUTE.zmin=BC(5);
FDTD.BoundaryCond.ATTRIBUTE.zmax=BC(6);
