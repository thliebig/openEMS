function FDTD = InitFDTD(NrTS, endCrit, varargin)

FDTD.ATTRIBUTE.NumberOfTimesteps=NrTS;
FDTD.ATTRIBUTE.endCriteria=endCrit;

for n=1:(nargin-2)/2
    FDTD.ATTRIBUTE.(varargin{2*n-1})=varargin{2*n};
end

