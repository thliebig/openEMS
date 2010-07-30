function FDTD = SetBoundaryCond(FDTD,BC)
% FDTD = SetBoundaryCond(FDTD,BC)
%
% BC = [xmin xmax ymin ymax zmin zmax];
% ?min/?max: 0=PEC 1=PMC 2=MUR-ABC 3=PML-ABC
% 
% openEMS matlab interface
% -----------------------
% author: Thorsten Liebig

if (numel(BC)~=6)
    error('openEMS:SetBoundaryCond','wrong number of boundary conditions');
end

if isnumeric(BC)
    FDTD.BoundaryCond.ATTRIBUTE.xmin=BC(1);
    FDTD.BoundaryCond.ATTRIBUTE.xmax=BC(2);
    FDTD.BoundaryCond.ATTRIBUTE.ymin=BC(3);
    FDTD.BoundaryCond.ATTRIBUTE.ymax=BC(4);
    FDTD.BoundaryCond.ATTRIBUTE.zmin=BC(5);
    FDTD.BoundaryCond.ATTRIBUTE.zmax=BC(6);
elseif iscell(BC)
    FDTD.BoundaryCond.ATTRIBUTE.xmin=BC{1};
    FDTD.BoundaryCond.ATTRIBUTE.xmax=BC{2};
    FDTD.BoundaryCond.ATTRIBUTE.ymin=BC{3};
    FDTD.BoundaryCond.ATTRIBUTE.ymax=BC{4};
    FDTD.BoundaryCond.ATTRIBUTE.zmin=BC{5};
    FDTD.BoundaryCond.ATTRIBUTE.zmax=BC{6};
else
    error('openEMS:SetBoundaryCond','unknown boundary condition type');
end

        