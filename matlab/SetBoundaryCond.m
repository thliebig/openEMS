function FDTD = SetBoundaryCond(FDTD, BC, varargin)
% FDTD = SetBoundaryCond(FDTD, BC, varargin)
%
% BC = [xmin xmax ymin ymax zmin zmax];
% or BC = {xmin xmax ymin ymax zmin zmax};
% ?min/?max: 
%   0 = PEC      or  'PEC'
%   1 = PMC      or  'PMC'
%   2 = MUR-ABC  or  'MUR'
%   3 = PML-ABC  or  'PML_x' with pml size x => 4..50
% 
% example:
% BC = [ 1     1     0     0     2     3     ]  %using numbers or
% BC = {'PMC' 'PMC' 'PEC' 'PEC' 'MUR' 'PML_8'}  %using equivalent strings
%
% mur-abc definitions
% define a phase-velocity to be used by the mur-abc
% useful e.g. for dispersive waveguides
% FDTD = SetBoundaryCond(FDTD,BC,'MUR_PhaseVelocity',299792457.93272);
% 
% 
% pml definitions
% 	arguments:  'PML_Grading','gradFunction'
% 		Define the pml grading grading function.
% 		Predefined variables in this grading function are:
% 			D  = depth in the pml in meter
% 			dl = mesh delta inside the pml in meter
% 			W  = width (length) of the pml in meter
% 			N  = number of cells for the pml
% 			Z  = wave impedance at the current depth and position
% 
% example: 
% FDTD = SetBoundaryCond(FDTD,BC); 
% or
% FDTD = SetBoundaryCond(FDTD,BC,'PML_Grading','-log(1e-6)*log(2.5)/(2*dl*pow(2.5,W/dl)-1) * pow(2.5, D/dl) / Z');
%
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


for n=1:(nargin-2)/2
    FDTD.BoundaryCond.ATTRIBUTE.(varargin{2*n-1}) = varargin{2*n};
end
        