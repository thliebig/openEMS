function FDTD = SetupMPI(FDTD, varargin)
% function FDTD = SetupMPI(FDTD, varargin);
%
% % example, split the FDTD mesh in 2 equal parts in x-direction
% % and split the FDTD mesh in 3 parts in z-direction, split at z=-500 and z=500
% % this will need a Settings.MPI.NrProc of 2*3=6
% FDTD = SetupMPI(FDTD,'SplitN_X',2 ,'SplitPos_Z', '-500,500');
% 
% See also RunOpenEMS_MPI
% 
% openEMS matlab interface
% -----------------------
% author: Thorsten Liebig

for n=1:(nargin-1)/2
    FDTD.MPI.ATTRIBUTE.(varargin{2*n-1})=varargin{2*n};
end
