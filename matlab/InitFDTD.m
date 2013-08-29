function FDTD = InitFDTD(varargin)
% function FDTD = InitFDTD(varargin)
%
% Inititalize the FDTD data-structure.
%
% optional field arguments for usage with openEMS:
%   NrTS:           max. number of timesteps to simulate (e.g. default=1e9)
%   EndCriteria:    end criteria, e.g. 1e-5, simulations stops if energy has
%                   decayed by this value (<1e-4 is recommended, default=1e-5)
%   MaxTime:        max. real time in seconds to simulate
%   OverSampling:   nyquist oversampling of time domain dumps
%   CoordSystem:    choose coordinate system (0 Cartesian, 1 Cylindrical)
%   MultiGrid:      define a cylindrical sub-grid radius
%   TimeStep:       force to use a given timestep (dangerous!)
%   TimeStepFactor: reduce the timestep by a given factor (>0 to <=1)
%   TimeStepMethod: 1 or 3 chose timestep method (1=CFL, 3=Rennigs (default))
%   CellConstantMaterial: set to 1 to assume a material is constant inside
%                         a cell (material probing in cell center)
%
% examples:
%   %default init with 1e9 max. timesteps and -50dB end-criteria
%   FDTD = InitFDTD();
%
%   %init with 1e6 max. timesteps and -60dB end-criteria
%   FDTD = InitFDTD('NrTS', 1e6, 'EndCriteria', 1e-6);
%
%   %cylindrical FDTD simulation
%   FDTD = InitFDTD('CoordSystem', 1);
%
%   See also InitCSX
%
% openEMS matlab interface
% -----------------------
% author: Thorsten Liebig (c) 2010-2013

% default values
NrTS = 1e9;
endCrit = 1e-5;

% legacy support
if ((nargin==1) && (isnumeric(varargin{1})))
    NrTS = varargin{1};
    warning('openEMS:InitFDTD',['Syntax for InitFDTD has changed, use: "InitFDTD(''NrTS'', ' num2str(NrTS) ')" instead! Legacy support enabled.']);
elseif ((nargin>1) && (isnumeric(varargin{1})) && (isnumeric(varargin{2})))
    NrTS = varargin{1};
    endCrit = varargin{2};
    varargin(1:2) = [];
    warning('openEMS:InitFDTD',['Syntax for InitFDTD has changed, use: "InitFDTD(''NrTS'', ' num2str(NrTS) ', ''EndCriteria'', ' num2str(endCrit) ')" instead! Legacy support enabled.']);
end

for n=1:numel(varargin)/2
    if strcmpi(varargin{2*n-1},'CoordSystem')==1
        FDTD.ATTRIBUTE.CylinderCoords=varargin{2*n}==1;
    elseif strcmpi(varargin{2*n-1},'NrTS')==1
        NrTS=varargin{2*n};
    elseif strcmpi(varargin{2*n-1},'EndCriteria')==1
        endCrit=varargin{2*n};
    else
        FDTD.ATTRIBUTE.(varargin{2*n-1})=varargin{2*n};
    end
end

FDTD.ATTRIBUTE.NumberOfTimesteps=NrTS;
FDTD.ATTRIBUTE.endCriteria=endCrit;
