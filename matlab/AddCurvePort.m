function [CSX,port] = AddCurvePort( CSX, prio, portnr, R, start, stop, excitename )
%[CSX,port] = AddCurvePort( CSX, prio, portnr, R, start, stop [, excitename] )
%
% Creates a curve port (1-dimensional).
% The mesh must already be initialized.
%
% input:
%   CSX:    CSX-object created by InitCSX()
%   prio:   priority for excitation, metal, sheet and probe boxes
%   portnr: (integer) number of the port
%   R:      internal resistance of the port
%   start:  3D start rowvector for port definition
%   stop:   3D end rowvector for port definition
%   excitename (optional): if specified, the port will be switched on (see AddExcitation())
%
% output:
%   CSX:
%   port:
%
% example:
%   start = [0 0 0]; stop = [0 0 12];
%   this defines a lumped port in z-direction
%   the excitation/probe is placed between start(1) and stop(1)
%
% (C) 2010 Sebastian Held <sebastian.held@uni-due.de>
% See also InitCSX AddExcitation

% make row vector
start = reshape( start, 1, [] );
stop  = reshape( stop , 1, [] );

% get grid
mesh{1} = sort(unique(CSX.RectilinearGrid.XLines));
mesh{2} = sort(unique(CSX.RectilinearGrid.YLines));
mesh{3} = sort(unique(CSX.RectilinearGrid.ZLines));
unit    = CSX.RectilinearGrid.ATTRIBUTE.DeltaUnit;

% find port direction
dir = abs(stop - start);
[dummy,dir] = max(dir);

% other directions
dir1 = mod(dir,3)+1;
dir2 = mod(dir+1,3)+1;

% normalize start and stop
if start(dir) < stop(dir)
    nstart = start;
    nstop  = stop;
else
    nstart = stop;
    nstop  = start;
end

% snap to grid
start_idx = zeros(1,3);
stop_idx  = zeros(1,3);
for n=1:3
    start_idx(n) = interp1( mesh{n}, 1:numel(mesh{n}), nstart(n), 'nearest' );
    stop_idx(n)  = interp1( mesh{n}, 1:numel(mesh{n}), nstop(n), 'nearest' );
end

% calculate position
port_start_idx = start_idx;
port_stop_idx  = stop_idx;
if abs(start_idx(dir) - stop_idx(dir)) ~= 1
    % calc port position
    idx  = interp1( mesh{dir},  1:numel(mesh{dir}),  (nstart(dir)+nstop(dir))/2,   'nearest' );
    idx1 = interp1( mesh{dir1}, 1:numel(mesh{dir1}), (nstart(dir1)+nstop(dir1))/2, 'nearest' );
    idx2 = interp1( mesh{dir2}, 1:numel(mesh{dir2}), (nstart(dir2)+nstop(dir2))/2, 'nearest' );
    port_start_idx(dir)  = idx;
    port_start_idx(dir1) = idx1;
    port_start_idx(dir2) = idx2;
    port_stop_idx(dir)   = idx+1;
    port_stop_idx(dir1)  = idx1;
    port_stop_idx(dir2)  = idx2;
    metalname = ['port' num2str(portnr) '_PEC'];
    CSX = AddMetal( CSX, metalname );
    CSX = AddCurve( CSX, metalname, prio, [nstart.' [mesh{1}(port_start_idx(1));mesh{2}(port_start_idx(2));mesh{3}(port_start_idx(3))]] );
    CSX = AddCurve( CSX, metalname, prio, [nstop.' [mesh{1}(port_stop_idx(1));mesh{2}(port_stop_idx(2));mesh{3}(port_stop_idx(3))]] );
end

% calculate position of resistive material
delta1_n = mesh{dir1}(port_start_idx(dir1)) - mesh{dir1}(port_start_idx(dir1)-1);
delta1_p = mesh{dir1}(port_start_idx(dir1)+1) - mesh{dir1}(port_start_idx(dir1));
delta2_n = mesh{dir2}(port_start_idx(dir2)) - mesh{dir2}(port_start_idx(dir2)-1);
delta2_p = mesh{dir2}(port_start_idx(dir2)+1) - mesh{dir2}(port_start_idx(dir2));
m_start = zeros(1,3);
m_stop  = zeros(1,3);
for n=1:3
    m_start(n)  = mesh{n}(port_start_idx(n));
    m_stop(n)   = mesh{n}(port_stop_idx(n));
end
m_start(dir1) = m_start(dir1) - delta1_n/2;
m_stop(dir1)  = m_stop(dir1)  + delta1_p/2;
m_start(dir2) = m_start(dir2) - delta2_n/2;
m_stop(dir2)  = m_stop(dir2)  + delta2_p/2;

% calculate kappa
l = abs((m_stop(dir) - m_start(dir))) * unit; % length of the "resistor"
A = (m_stop(dir1) - m_start(dir1)) * (m_stop(dir2) - m_start(dir2)) * unit^2; % area of the "resistor"
A = abs(A);
kappa = l/A / R; % [kappa] = S/m
materialname = ['port' num2str(portnr) '_sheet_resistance'];
CSX = AddMaterial( CSX, materialname );%, 'Isotropy', 0 );
% kappa_cell = {0,0,0};
% kappa_cell{dir} = kappa;
kappa_cell = kappa;
CSX = SetMaterialProperty( CSX, materialname, 'Kappa', kappa_cell );
CSX = AddBox( CSX, materialname, prio, m_start, m_stop );

% calculate position of the voltage probe
v_start = [mesh{1}(port_start_idx(1)), mesh{2}(port_start_idx(2)), mesh{3}(port_start_idx(3))];
v_stop  = [mesh{1}(port_stop_idx(1)),  mesh{2}(port_stop_idx(2)),  mesh{3}(port_stop_idx(3))];

% calculate position of the current probe
i_start = m_start;
i_stop  = m_stop;
i_start(dir) = (i_start(dir)+i_stop(dir))/2;
i_stop(dir)  = i_start(dir);

% create the probes
name = ['port_ut' num2str(portnr)];
weight = -1;
CSX = AddProbe( CSX, name, 0, weight );
CSX = AddBox( CSX, name, prio, v_start, v_stop );
name = ['port_it' num2str(portnr)];
weight = 1;
CSX = AddProbe( CSX, name, 1, weight );
CSX = AddBox( CSX, name, prio, i_start, i_stop );

% create port structure
port.nr = portnr;
port.drawingunit = unit;
% port.start = start;
% port.stop = stop;
% port.v_start = v_start;
% port.v_stop  = v_stop;
% port.i_start = i_start;
% port.i_stop  = i_stop;
% port.dir = dir;
% port.direction = direction;
% port.idx_cal = idx_cal;
% port.idx1 = idx1;
% port.idx1 = idx1;
port.excite = 0;

% create excitation
if (nargin >= 6) && ~isempty(excitename)
	% excitation of this port is enabled
	port.excite = 1;
    e_start = v_start;
    e_stop  = v_stop;
    CSX = AddExcitation( CSX, excitename, 0, start_idx ~= stop_idx );
	CSX = AddBox( CSX, excitename, prio, e_start, e_stop );
end
