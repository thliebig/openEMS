function FDTD = SetupMPI(FDTD, varargin)
% function FDTD = SetupMPI(FDTD, varargin);
%
% % example
% FDTD = SetupMPI(FDTD,'SplitPos_X', '-500,500','SplitPos_Z', '-500,500');
% 
% openEMS matlab interface
% -----------------------
% author: Thorsten Liebig

for n=1:(nargin-1)/2
    FDTD.MPI.ATTRIBUTE.(varargin{2*n-1})=varargin{2*n};
end
