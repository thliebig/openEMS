function [CSX,port] = AddWaveGuidePort( CSX, prio, portnr, start, stop, dir, E_WG_func, H_WG_func, kc, exc_amp, varargin )
% function [CSX,port] = AddWaveGuidePort( CSX, prio, portnr, start, stop, dir, E_WG_func, H_WG_func, kc, exc_amp, varargin )
% 
% Create a waveguide port, including an optional excitation and probes
% 
% Note: - The excitation will be located at the start position in the given direction
%       - The voltage and current probes at the stop position in the given direction
%
% parameter:
%   CSX:        complete CSX structure (must contain a mesh)
%   prio:       priority of primitives
%   start:      start coordinates of waveguide port box
%   stop:       stop  coordinates of waveguide port box
%   dir:        direction of port (0/1/2 or 'x'/'y'/'z'-direction)
%   E_WG_func:  electric field mode profile function as a string
%   H_WG_func:  magnetic field mode profile function as a string
%   kc:         cutoff wavenumber (defined by the waveguide dimensions)
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
%   p11 = 1.841;
%   kc = p11 / radius;  % cutoff wavenumber with radius in meter
%   kc_draw = kc*unit;  % cutoff wavenumber in drawing units
%
%   % electric field mode profile
%   func_E{1} = [ num2str(-1/kc_draw^2,15) '/rho*cos(a)*j1('  num2str(kc_draw,15) '*rho)'];
%   func_E{2} = [ num2str(1/kc_draw,15) '*sin(a)*0.5*(j0('  num2str(kc_draw,15) '*rho)-jn(2,'  num2str(kc_draw,15) '*rho))'];
%   func_E{3} = 0;
%
%   % magnetic field mode profile
%   func_H{1} = [ '-1*' num2str(1/kc_draw,15) '*sin(a)*0.5*(j0('  num2str(kc_draw,15) '*rho)-jn(2,'  num2str(kc_draw,15) '*rho))'];
%   func_H{2} = [ num2str(-1/kc_draw^2,15) '/rho*cos(a)*j1('  num2str(kc_draw,15) '*rho)'];
%   func_H{3} = 0;
%
%   start=[mesh.r(1)   mesh.a(1)   0  ];
%   stop =[mesh.r(end) mesh.a(end) 100];
%   [CSX, port{1}] = AddWaveGuidePort(CSX, 0, 1, start, stop, 2, func_E, func_H, kc, 1);
%
% openEMS matlab interface
% -----------------------
% (c) 2013 Thorsten Liebig (thorsten.liebig@gmx.de)
%
% See also InitCSX, AddExcitation, calcWGPort, calcPort

%check mesh
if ~isfield(CSX,'RectilinearGrid')
    error 'mesh needs to be defined! Use DefineRectGrid() first!';
end

dir = DirChar2Int(dir);

port.type='WaveGuide';
port.nr=portnr;
port.kc = kc;
port.dir = dir;
port.drawingunit = CSX.RectilinearGrid.ATTRIBUTE.DeltaUnit;

PortNamePrefix = '';

varargin_tmp  = varargin;
for n=1:2:numel(varargin_tmp)
    if strcmpi('PortNamePrefix',varargin_tmp{n})
        PortNamePrefix = varargin_tmp{n+1};
        varargin([n n+1]) = [];
    end
end

% matlab adressing
dir = dir + 1;
dir_sign = sign(stop(dir) - start(dir));
if (dir_sign==0)
    dir_sign = 1;
end

port.direction = dir_sign;

E_WG_func{dir} = 0;
H_WG_func{dir} = 0;

port.excite = 0;
if (exc_amp~=0)
    if (start(dir)==stop(dir))
        error 'if waveguide port is to be excited, the length in propagation direction must not be zero'
    end
    e_start = start;
    e_stop = stop;
    e_stop(dir) = e_start(dir);
    port.excite = 1;
    port.excitepos = e_start(dir);
    e_vec = [1 1 1]*exc_amp;
    e_vec(dir) = 0;
    exc_name = [PortNamePrefix 'port_excite_' num2str(portnr)];
    CSX = AddExcitation( CSX, exc_name, 0, e_vec, varargin{:});
    CSX = SetExcitationWeight(CSX, exc_name, E_WG_func );
	CSX = AddBox( CSX, exc_name, prio, e_start, e_stop);
end

% voltage/current planes
m_start = start;
m_stop = stop;
m_start(dir) = stop(dir);

port.measplanepos = m_start(dir);
port.U_filename = [PortNamePrefix 'port_ut' int2str(portnr)];
CSX = AddProbe(CSX, port.U_filename, 10, 'ModeFunction', E_WG_func);
CSX = AddBox(CSX, port.U_filename, 0 ,m_start, m_stop);

port.I_filename = [PortNamePrefix 'port_it' int2str(portnr)];
CSX = AddProbe(CSX, port.I_filename, 11, 'ModeFunction', H_WG_func, 'weight', dir_sign);
CSX = AddBox(CSX, port.I_filename, 0 ,m_start, m_stop);
