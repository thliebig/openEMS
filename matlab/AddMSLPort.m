function [CSX,port] = AddMSLPort( CSX, portnr, materialname, start, stop, dir, evec, refplaneshift, excitename )
% [CSX,port] = AddMSLPort( CSX, portnr, materialname, start, stop, dir, evec, refplaneshift, excitename )
%
% CSX: CSX-object created by InitCSX()
% portnr: (integer) number of the port
% materialname: property for the MSL (created by AddMetal() or AddMaterial())
% start: 3D start rowvector for port definition
% stop:  3D end rowvector for port definition
% dir: direction of wave propagation (choices: [1 0 0], [0 1 0] or [0 0 1])
% evec: excitation vector, which defines the direction of the e-field (must be the same as used in AddExcitation())
% refplaneshift (optional): if not specified or empty, the measurement
%    plane is used; if specified, reference plane is shifted by
%    <refplaneshift> starting from <start> (thus refplaneshift is normally
%    positive)
% excitename (optional): if specified, the port will be switched on (see AddExcitation())
%
% the mesh must be already initialized
%
% example:
%   start = [0 0 height]; stop = [length width 0]; dir = [1 0 0]; evec = [0 0 -1]
%   this defines a MSL in x-direction (dir) with an e-field excitation in -z-direction (evec)
%   the excitation is placed at x=start(1); the wave travels towards x=stop(1)
%   the MSL-metal is created in xy-plane at z=start(3)
%
% Sebastian Held <sebastian.held@gmx.de>
% May 13 2010
%
% See also InitCSX AddMetal AddMaterial AddExcitation calcMSLPort

% check dir
if ~(dir(1) == dir(2) == 0) && ~(dir(1) == dir(3) == 0) && ~(dir(2) == dir(3) == 0) || (sum(dir) == 0)
	error 'dir must have exactly one component ~= 0'
end
dir = dir ./ sum(dir); % dir is now a unit vector

% check evec
if ~(evec(1) == evec(2) == 0) && ~(evec(1) == evec(3) == 0) && ~(evec(2) == evec(3) == 0) || (sum(evec) == 0)
	error 'evec must have exactly one component ~= 0'
end
evec0 = evec ./ sum(evec); % evec0 is a unit vector

% normalize start and stop
nstart = min( [start;stop] );
nstop  = max( [start;stop] );

% determine index (1, 2 or 3) of propagation (length of MSL)
idx_prop = dir * [1;2;3];

% determine index (1, 2 or 3) of width of MSL
idx_width = abs(cross(dir,evec0)) * [1;2;3];

% determine index (1, 2 or 3) of height
idx_height = abs(evec0) * [1;2;3];

% direction of propagation
if stop(idx_prop)-start(idx_prop) > 0
    direction = +1;
else
    direction = -1;
end

% create the metal/material for the MSL
MSL_start = start;
MSL_stop = stop;
MSL_stop(idx_height) = MSL_start(idx_height);
CSX = AddBox( CSX, materialname, 999, MSL_start, MSL_stop );

% FIXME
% openEMS v0.0.7 does not snap PEC

% calculate position of the voltage probes
mesh{1} = sort(CSX.RectilinearGrid.XLines);
mesh{2} = sort(CSX.RectilinearGrid.YLines);
mesh{3} = sort(CSX.RectilinearGrid.ZLines);
meshlines = interp1( mesh{idx_prop}, 1:numel(mesh{idx_prop}), (nstart(idx_prop)+nstop(idx_prop))/2, 'nearest' );
meshlines = mesh{idx_prop}(meshlines-1:meshlines+1); % get three lines (approx. at center)
if direction == -1
    meshlines = fliplr(meshlines);
end
MSL_w2 = interp1( mesh{idx_width}, 1:numel(mesh{idx_width}), (nstart(idx_width)+nstop(idx_width))/2, 'nearest' );
MSL_w2 = mesh{idx_width}(MSL_w2); % get e-line at center of MSL (MSL_width/2)
v1_start(idx_prop)   = meshlines(1);
v1_start(idx_width)  = MSL_w2;
v1_start(idx_height) = nstop(idx_height);
v1_stop  = v1_start;
v1_stop(idx_height)  = nstart(idx_height);
v2_start = v1_start;
v2_stop  = v1_stop;
v2_start(idx_prop)   = meshlines(2);
v2_stop(idx_prop)    = meshlines(2);
v3_start = v2_start;
v3_stop  = v2_stop;
v3_start(idx_prop)   = meshlines(3);
v3_stop(idx_prop)    = meshlines(3);

% calculate position of the current probes
idx = interp1( mesh{idx_width}, 1:numel(mesh{idx_width}), nstart(idx_width), 'nearest' );
i1_start(idx_width)  = mesh{idx_width}(idx) - diff(mesh{idx_width}(idx-1:idx))/2;
idx = interp1( mesh{idx_height}, 1:numel(mesh{idx_height}), start(idx_height), 'nearest' );
i1_start(idx_height) = mesh{idx_height}(idx) - diff(mesh{idx_height}(idx-1:idx))/2;
i1_stop(idx_height)  = mesh{idx_height}(idx) + diff(mesh{idx_height}(idx:idx+1))/2;
i1_start(idx_prop)   = sum(meshlines(1:2))/2;
i1_stop(idx_prop)    = i1_start(idx_prop);
idx = interp1( mesh{idx_width}, 1:numel(mesh{idx_width}), nstop(idx_width), 'nearest' );
i1_stop(idx_width)   = mesh{idx_width}(idx) + diff(mesh{idx_width}(idx:idx+1))/2;
i2_start = i1_start;
i2_stop  = i1_stop;
i2_start(idx_prop)   = sum(meshlines(2:3))/2;
i2_stop(idx_prop)    = i2_start(idx_prop);

% create the probes
name = ['port_ut' num2str(portnr) 'A'];
weight = sum(evec);
CSX = AddProbe( CSX, name, 0, weight );
CSX = AddBox( CSX, name, 999, v1_start, v1_stop );
name = ['port_ut' num2str(portnr) 'B'];
CSX = AddProbe( CSX, name, 0, weight );
CSX = AddBox( CSX, name, 999, v2_start, v2_stop );
name = ['port_ut' num2str(portnr) 'C'];
CSX = AddProbe( CSX, name, 0, weight );
CSX = AddBox( CSX, name, 999, v3_start, v3_stop );
name = ['port_it' num2str(portnr) 'A'];
weight = direction;
CSX = AddProbe( CSX, name, 1, weight );
CSX = AddBox( CSX, name, 999, i1_start, i1_stop );
name = ['port_it' num2str(portnr) 'B'];
CSX = AddProbe( CSX, name, 1, weight );
CSX = AddBox( CSX, name, 999, i2_start, i2_stop );

% create port structure
port.nr = portnr;
port.drawingunit = CSX.RectilinearGrid.ATTRIBUTE.DeltaUnit;
% port.start = start;
% port.stop = stop;
% port.v1_start = v1_start;
% port.v1_stop = v1_stop;
% port.v2_start = v2_start;
% port.v2_stop = v2_stop;
% port.v3_start = v3_start;
% port.v3_stop = v3_stop;
port.v_delta = diff(meshlines);
% port.i1_start = i1_start;
% port.i1_stop = i1_stop;
% port.i2_start = i2_start;
% port.i2_stop = i2_stop;
port.i_delta = diff( meshlines(1:end-1) + diff(meshlines)/2 );
port.direction = direction;
% port.dir = dir;
% port.evec = evec;
% port.idx_prop = idx_prop;
% port.idx_width = idx_width;
% port.idx_height = idx_height;
port.excite = 0;
port.refplaneshift = 0;
port.measplanepos = abs(v2_start(idx_prop) - start(idx_prop));

if (nargin >= 8) && (~isempty(refplaneshift))
    % refplaneshift counts from start of port
    port.refplaneshift = refplaneshift - direction*(v2_start(idx_prop) - start(idx_prop));
end

% create excitation
if nargin >= 9
	% excitation of this port is enabled
	port.excite = 1;
    meshline = interp1( mesh{idx_prop}, 1:numel(mesh{idx_prop}), start(idx_prop), 'nearest' );
    ex_start(idx_prop)   = mesh{idx_prop}(meshline+direction);
	ex_start(idx_width)  = nstart(idx_width);
	ex_start(idx_height) = nstart(idx_height);
	ex_stop(idx_prop)    = ex_start(idx_prop);
	ex_stop(idx_width)   = nstop(idx_width);
	ex_stop(idx_height)  = nstop(idx_height);
    CSX = AddExcitation( CSX, excitename, 0, evec );
	CSX = AddBox( CSX, excitename, 999, ex_start, ex_stop );
end
