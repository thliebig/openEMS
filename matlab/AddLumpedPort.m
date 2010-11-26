function [CSX,port] = AddLumpedPort( CSX, portnr, R, start, stop, dir, excitename )
% [CSX,port] = AddLumpedPort( CSX, portnr, materialname, start, stop, dir, evec, excitename )
%
% CSX: CSX-object created by InitCSX()
% portnr: (integer) number of the port
% R: internal resistance of the port
% start: 3D start rowvector for port definition
% stop:  3D end rowvector for port definition
% dir: direction of of port (choices: [1 0 0], [0 1 0] or [0 0 1])
% excitename (optional): if specified, the port will be switched on (see AddExcitation())
%
% the mesh must be already initialized
%
% example:
%   start = [0 0 height]; stop = [length width height]; dir = [1 0 0];
%   this defines a lumped port in x-direction (dir)
%   the excitation/probe is placed between start(1) and stop(1)
%
% Sebastian Held <sebastian.held@gmx.de>
% Jun 1 2010
%
% See also InitCSX AddExcitation

% check dir
if ~(dir(1) == dir(2) == 0) && ~(dir(1) == dir(3) == 0) && ~(dir(2) == dir(3) == 0) || (sum(dir) == 0)
	error 'dir must have exactly one component ~= 0'
end
dir = dir ./ sum(dir); % dir is now a unit vector

% get grid
mesh{1} = sort(CSX.RectilinearGrid.XLines);
mesh{2} = sort(CSX.RectilinearGrid.YLines);
mesh{3} = sort(CSX.RectilinearGrid.ZLines);
drawingunit = CSX.RectilinearGrid.ATTRIBUTE.DeltaUnit;

% snap to grid
idx_plane = 0;
for n=1:3
    start_idx = interp1( mesh{n}, 1:numel(mesh{n}), start(n), 'nearest' );
    stop_idx  = interp1( mesh{n}, 1:numel(mesh{n}), stop(n),  'nearest' );
    if start_idx == stop_idx
        idx_plane = n; % two dimensional port: this is the correct index
    end
    start(n)  = mesh{n}(start_idx);
    stop(n)   = mesh{n}(stop_idx);
end

if idx_plane == 0
    error( 'the port must be two-dimensional!' );
end

% normalize start and stop
nstart = min( [start;stop] );
nstop  = max( [start;stop] );

% determine index (1, 2 or 3) of calibration (e-) line
idx_cal = dir * [1;2;3];

% direction of calibration line
if stop(idx_cal)-start(idx_cal) > 0
    direction = +1;
else
    direction = -1;
end

% determine the other direction (FIXME is there a better way?)
idx1 = [1 2 3];
idx1 = idx1(idx1 ~= idx_plane);
idx1 = idx1(idx1 ~= idx_cal);

% calculate position of resistive material
idx = interp1( mesh{idx_plane}, 1:numel(mesh{idx_plane}), nstart(idx_plane), 'nearest' );
delta2_n = mesh{idx_plane}(idx) - mesh{idx_plane}(idx-1);
delta2_p = mesh{idx_plane}(idx+1) - mesh{idx_plane}(idx);
m_start = nstart;
m_stop  = nstop;
m_start(idx_plane) = m_start(idx_plane) - delta2_n/2;
m_stop(idx_plane)  = m_stop(idx_plane)  + delta2_p/2;

% calculate kappa
l = (m_stop(idx_cal) - m_start(idx_cal)) * drawingunit; % length of the "sheet"
A = (m_stop(idx1) - m_start(idx1)) * (m_stop(idx_plane) - m_start(idx_plane)) * drawingunit^2; % area of the "sheet"
kappa = l/A / R; % [kappa] = S/m
CSX = AddMaterial( CSX, ['port' num2str(portnr) '_sheet_resistance'], 'Isotropy', 0 );
kappa_cell = {};
kappa_cell{1} = kappa*dir(1);
kappa_cell{2} = kappa*dir(2);
kappa_cell{3} = kappa*dir(3);
CSX = SetMaterialProperty( CSX, ['port' num2str(portnr) '_sheet_resistance'], 'Kappa', kappa_cell );
CSX = AddBox( CSX, ['port' num2str(portnr) '_sheet_resistance'], 999, m_start, m_stop );

% calculate position of the voltage probe
v_start(idx_plane) = start(idx_plane);
center1 = interp1( mesh{idx1}, 1:numel(mesh{idx1}), (nstart(idx1)+nstop(idx1))/2, 'nearest' );
v_start(idx1) = mesh{idx1}(center1);
v_stop = v_start;
v_start(idx_cal) = nstart(idx_cal);
v_stop(idx_cal)  = nstop(idx_cal);

% calculate position of the current probe
idx = interp1( mesh{idx1}, 1:numel(mesh{idx1}), nstart(idx1), 'nearest' );
delta1_n = mesh{idx1}(idx) - mesh{idx1}(idx-1);
idx = interp1( mesh{idx1}, 1:numel(mesh{idx1}), nstop(idx1), 'nearest' );
delta1_p = mesh{idx1}(idx+1) - mesh{idx1}(idx);
h_offset = diff(mesh{idx_cal});
idx = interp1( mesh{idx_cal} + [h_offset h_offset(end)]/2, 1:numel(mesh{idx_cal}), (nstart(idx_cal)+nstop(idx_cal))/2, 'nearest' );
i_start(idx_cal)   = mesh{idx_cal}(idx) + h_offset(idx)/2;
i_stop(idx_cal)    = i_start(idx_cal);
i_start(idx1)      = nstart(idx1) - delta1_n/2;
i_start(idx_plane) = nstart(idx_plane) - delta2_n/2;
i_stop(idx1)       = nstop(idx1) + delta1_p/2;
i_stop(idx_plane)  = nstop(idx_plane) + delta2_p/2;

% create the probes
name = ['port_ut' num2str(portnr)];
weight = -direction;
CSX = AddProbe( CSX, name, 0, weight );
CSX = AddBox( CSX, name, 999, v_start, v_stop );
name = ['port_it' num2str(portnr)];
weight = direction;
CSX = AddProbe( CSX, name, 1, weight );
CSX = AddBox( CSX, name, 999, i_start, i_stop );

% create port structure
port.nr = portnr;
port.drawingunit = CSX.RectilinearGrid.ATTRIBUTE.DeltaUnit;
% port.start = start;
% port.stop = stop;
% port.v_start = v_start;
% port.v_stop  = v_stop;
% port.i_start = i_start;
% port.i_stop  = i_stop;
% port.dir = dir;
port.direction = direction;
% port.idx_cal = idx_cal;
% port.idx1 = idx1;
% port.idx1 = idx1;
port.excite = 0;

% create excitation
if (nargin >= 7) && ~isempty(excitename)
	% excitation of this port is enabled
	port.excite = 1;
    e_start = nstart;
    e_stop  = nstop;
    e_start(idx_plane) = start(idx_plane); % excitation-plane is determined by start vector
    e_stop(idx_plane)  = start(idx_plane);
    CSX = AddExcitation( CSX, excitename, 0, -dir*direction);
	CSX = AddBox( CSX, excitename, 999, e_start, e_stop );
end
