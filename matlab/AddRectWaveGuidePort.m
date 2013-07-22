function [CSX,port] = AddRectWaveGuidePort( CSX, prio, portnr, start, stop, dir, a, b, mode_name, exc_amp, varargin )
% function [CSX,port] = AddRectWaveGuidePort( CSX, prio, portnr, start, stop, dir, a, b, mode_name, exc_amp, varargin )
% 
% Create a rectangular waveguide port, including an optional excitation and probes
% 
% Note: - The excitation will be located at the start position in the given direction
%       - The voltage and current probes at the stop position in the given direction
%
% input:
%   CSX:        complete CSX structure (must contain a mesh)
%   prio:       priority of primitives
%   start:      start coordinates of waveguide port box
%   stop:       stop  coordinates of waveguide port box
%   dir:        direction of port (0/1/2 or 'x'/'y'/'z'-direction)
%   a,b:        rectangular waveguide width and height (in meter)
%   mode_name:  mode name, e.g. 'TE11' or 'TM21'
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
%   % create a TE10 circular waveguide mode, using cylindircal coordinates
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

if (nargin<10)
    exc_amp = 0;
end

m = str2double(mode_name(3));
n = str2double(mode_name(4));

% values by David M. Pozar, Microwave Engineering, third edition
kc = sqrt((m*pi/a)^2 + (n*pi/b)^2);

if ~isfield(CSX,'RectilinearGrid')
    error 'mesh needs to be defined! Use DefineRectGrid() first!';
    if (~isfield(CSX.RectilinearGrid,'XLines') || ~isfield(CSX.RectilinearGrid,'YLines') || ~isfield(CSX.RectilinearGrid,'ZLines'))
        error 'mesh needs to be defined! Use DefineRectGrid() first!';
    end
end

unit = CSX.RectilinearGrid.ATTRIBUTE.DeltaUnit;

dir = DirChar2Int(dir);
dir_names={'x','y','z'};

dirP = mod((dir+1),3)+1;
dirPP = mod((dir+2),3)+1;
nameX = ['(' dir_names{dirP}  '-' num2str(start(dirP)) ')'];
nameY = ['(' dir_names{dirPP} '-' num2str(start(dirPP)) ')'];

%convert a&b to drawing units
a = a/unit;
b = b/unit;
% functions by David M. Pozar, Microwave Engineering, third edition
% electric field mode profile
func_Ex = [num2str( n/b) '*cos(' num2str(m*pi/a) '*' nameX ')*sin('  num2str(n*pi/b) '*' nameY ')'];
func_Ey = [num2str(-m/a) '*sin(' num2str(m*pi/a) '*' nameX ')*cos('  num2str(n*pi/b) '*' nameY ')'];

% magnetic field mode profile
func_Hx = [num2str(m/a) '*sin(' num2str(m*pi/a) '*' nameX ')*cos('  num2str(n*pi/b) '*' nameY ')'];
func_Hy = [num2str(n/b) '*cos(' num2str(m*pi/a) '*' nameX ')*sin('  num2str(n*pi/b) '*' nameY ')'];


func_E{dir+1} = 0;
func_E{dirP} = func_Ex;
func_E{dirPP} = func_Ey;

func_H{dir+1} = 0;
func_H{dirP} = func_Hx;
func_H{dirPP} = func_Hy;

[CSX,port] = AddWaveGuidePort( CSX, prio, portnr, start, stop, dir, func_E, func_H, kc, exc_amp, varargin{:} );

