function [CSX,port] = AddCurvePort( CSX, prio, portnr, R, start, stop, excite, varargin )
%[CSX,port] = AddCurvePort( CSX, prio, portnr, R, start, stop [, excite, varargin] )
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
%   excite (optional): if true, the port will be switched on (see AddExcitation())
%                       Note: for legacy support a string will be accepted
% optional (key/values):
%   varargin:   optional additional excitations options, see also AddExcitation
%   'PortNamePrefix': a prefix to the port name
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

port.type='Lumped';
port.nr=portnr;

PortNamePrefix = '';

varargin_tmp  = varargin;
for n=1:2:numel(varargin_tmp)
    if strcmpi('PortNamePrefix',varargin_tmp{n})
        PortNamePrefix = varargin_tmp{n+1};
        varargin([n n+1]) = [];
    end
end

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
try
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
		metalname = [PortNamePrefix 'port' num2str(portnr) '_PEC'];
		CSX = AddMetal( CSX, metalname );
		CSX = AddCurve( CSX, metalname, prio, [nstart.' [mesh{1}(port_start_idx(1));mesh{2}(port_start_idx(2));mesh{3}(port_start_idx(3))]] );
		CSX = AddCurve( CSX, metalname, prio, [nstop.' [mesh{1}(port_stop_idx(1));mesh{2}(port_stop_idx(2));mesh{3}(port_stop_idx(3))]] );
	end
catch
	error('Unable to place port on mesh; check the location of the port, and make sure that the mesh is large enough');
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

% calculate position of the voltage probe & excitation
v_start = [mesh{1}(port_start_idx(1)), mesh{2}(port_start_idx(2)), mesh{3}(port_start_idx(3))];
v_stop  = [mesh{1}(port_stop_idx(1)),  mesh{2}(port_stop_idx(2)),  mesh{3}(port_stop_idx(3))];

% calculate position of the current probe
i_start = m_start;
i_stop  = m_stop;
i_start(dir) = (i_start(dir)+i_stop(dir))/2;
i_stop(dir)  = i_start(dir);

% create the probes
port.U_filename = [PortNamePrefix 'port_ut' num2str(portnr)];
weight = -1;
CSX = AddProbe( CSX, port.U_filename, 0, 'weight', weight );
CSX = AddBox( CSX, port.U_filename, prio, v_start, v_stop );
port.I_filename = [PortNamePrefix 'port_it' num2str(portnr)];
weight = 1;
CSX = AddProbe( CSX, port.I_filename, 1, 'weight',  weight );
CSX = AddBox( CSX, port.I_filename, prio, i_start, i_stop );

% create port structure
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

if (nargin < 7)
    excite = false;
end

% legacy support, will be removed at some point
if ischar(excite)
    warning('CSXCAD:AddCurvePort','depreceated: a string as excite option is no longer supported and will be removed in the future, please use true or false');
    if ~isempty(excite)
        excite = true;
    else
        excite = false;
    end
end

port.excite = excite;

% create excitation
if (excite)
	% excitation of this port is enabled
	port.excite = 1;
    e_start = v_start;
    e_stop  = v_stop;
    CSX = AddExcitation( CSX, [PortNamePrefix 'port_excite_' num2str(portnr)], 0, start_idx ~= stop_idx, varargin{:});
	CSX = AddBox( CSX, [PortNamePrefix 'port_excite_' num2str(portnr)], prio, e_start, e_stop );
end

port.Feed_R = R;
if (R>0 && (~isinf(R))) 
    CSX = AddLumpedElement( CSX, [PortNamePrefix 'port_resist_' int2str(portnr)], dir-1, 'R', R);
    CSX = AddBox( CSX, [PortNamePrefix 'port_resist_' int2str(portnr)], prio, v_start, v_stop );
elseif (R==0)
    CSX = AddBox(CSX,metalname, prio, v_start, v_stop);
end
