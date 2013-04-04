function [CSX,port] = AddCircWaveGuidePort( CSX, prio, portnr, start, stop, radius, mode_name, pol_ang, exc_amp, varargin )
% function [CSX,port] = AddCircWaveGuidePort( CSX, prio, portnr, start, stop, radius, mode_name, pol_ang, exc_amp, varargin )
% 
% Create a circular waveguide port, including an optional excitation and probes
% 
% Note: - The excitation will be located at the start position in the given direction
%       - The voltage and current probes at the stop position in the given direction
%
% input:
%   CSX:        complete CSX structure (must contain a mesh)
%   prio:       priority of primitives
%   start:      start coordinates of waveguide port box
%   stop:       stop  coordinates of waveguide port box
%   radius:     circular waveguide radius (in meter)
%   mode_name:  mode name, e.g. 'TE11' or 'TM21'
%   pol_ang:    polarization angle (e.g. 0 = horizontal, pi/2 = vertical)
%   exc_amp:    excitation amplitude (set 0 to be passive)
%
% optional (key/values):
%   varargin:   optional additional excitations options, see also AddExcitation
%   'PortNamePrefix': a prefix to the port name
%
% output:
%   CSX:        modified CSX structure
%   port:       port structure to use with calcPort
%
% example:
%   % create a TE11 circular waveguide mode, using cylindircal coordinates
%   start=[mesh.r(1)   mesh.a(1)   0  ];
%   stop =[mesh.r(end) mesh.a(end) 100];
%   [CSX,port] = AddCircWaveGuidePort( CSX, 99, 1, start, stop, 320e-3, 'TE11', 0, 1);
%
% openEMS matlab interface
% -----------------------
% (c) 2013 Thorsten Liebig (thorsten.liebig@gmx.de)
%
% See also InitCSX, AddExcitation, calcWGPort, calcPort

if (~strcmpi(mode_name(1:2),'TE'))
    error 'currently only TE type modes are supported'
end

if (nargin<9)
    exc_amp = 0;
end
if (nargin<8)
    pol_ang = 0;
end

pnm = 0;
n = str2double(mode_name(3));
m = str2double(mode_name(4));

% values by David M. Pozar, Microwave Engineering, third edition
if ((n==0) && (m==1))
    pnm = 3.832;
elseif ((n==1) && (m==1))
    pnm = 1.841;
elseif ((n==2) && (m==1))
    pnm = 3.054;
elseif ((n==0) && (m==2))
    pnm = 7.016;
elseif ((n==1) && (m==2))
    pnm = 5.331;
elseif ((n==2) && (m==2))
    pnm = 6.706;
elseif ((n==0) && (m==3))
    pnm = 10.174;
elseif ((n==1) && (m==3))
    pnm = 8.536;
elseif ((n==2) && (m==3))
    pnm = 9.970;
else
    error 'invalid TE_nm mode'
end

if ~isfield(CSX,'RectilinearGrid')
    error 'mesh needs to be defined! Use DefineRectGrid() first!';
end

unit = CSX.RectilinearGrid.ATTRIBUTE.DeltaUnit;
kc = pnm/radius;
kc_draw = kc*unit;

angle = ['a-' num2str(pol_ang)];
% functions by David M. Pozar, Microwave Engineering, third edition
% electric field mode profile
func_Er = [ num2str(-1/kc_draw^2,15) '/rho*cos(' angle ')*j1('  num2str(kc_draw,15) '*rho)'];
func_Ea = [ num2str(1/kc_draw,15) '*sin(' angle ')*0.5*(j0('  num2str(kc_draw,15) '*rho)-jn(2,'  num2str(kc_draw,15) '*rho))'];

% magnetic field mode profile
func_Hr = [ num2str(-1/kc_draw,15) '*sin(' angle ')*0.5*(j0('  num2str(kc_draw,15) '*rho)-jn(2,'  num2str(kc_draw,15) '*rho))'];
func_Ha = [ num2str(-1/kc_draw^2,15) '/rho*cos(' angle ')*j1('  num2str(kc_draw,15) '*rho)'];

if (CSX.ATTRIBUTE.CoordSystem==1)
    func_E = {func_Er, func_Ea, 0};
    func_H = {func_Hr, func_Ha, 0};
else
    func_Ex = ['(' func_Er '*cos(a) - ' func_Ea '*sin(a) ) * (rho<' num2str(radius/unit) ')'];
    func_Ey = ['(' func_Er '*sin(a) + ' func_Ea '*cos(a) ) * (rho<' num2str(radius/unit) ')'];
    func_E = {func_Ex, func_Ey, 0};
    
    func_Hx = ['(' func_Hr '*cos(a) - ' func_Ha '*sin(a) ) * (rho<' num2str(radius/unit) ')'];
    func_Hy = ['(' func_Hr '*sin(a) + ' func_Ha '*cos(a) ) * (rho<' num2str(radius/unit) ')'];
    func_H = {func_Hx, func_Hy, 0};
end

[CSX,port] = AddWaveGuidePort( CSX, prio, portnr, start, stop, 2, func_E, func_H, kc, exc_amp, varargin{:} );

