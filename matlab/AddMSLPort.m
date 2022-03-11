function [CSX,port] = AddMSLPort( CSX, prio, portnr, materialname, start, stop, dir, evec, varargin )
% [CSX,port] = AddMSLPort( CSX, prio, portnr, materialname, start, stop, dir, evec, varargin )
%
% CSX:          CSX-object created by InitCSX()
% prio:         priority for excitation and probe boxes
% portnr:       (integer) number of the port
% materialname: property for the MSL (created by AddMetal())
% start:        3D start rowvector for port definition
% stop:         3D end rowvector for port definition
% dir:          direction of wave propagation (choices: 0, 1, 2 or 'x','y','z')
% evec:         excitation vector, which defines the direction of the e-field (must be the same as used in AddExcitation())
% 
% variable input:
%  varargin:    optional additional excitations options, see also AddExcitation
% 'ExcitePort'  true/false to make the port an active feeding port (default
%               is false)
% 'FeedShift'   shift to port from start by a given distance in drawing
%               units. Default is 0. Only active if 'ExcitePort' is set!
% 'Feed_R'      Specify a lumped port resistance. Default is no lumped
%               port resistance --> port has to end in an ABC.
% 'MeasPlaneShift'  Shift the measurement plane from start t a given distance 
%               in drawing units. Default is the middle of start/stop.
% 'PortNamePrefix' a prefix to the port name
% 
% Important: The mesh has to be already set and defined by DefineRectGrid!
%
% example:
%   CSX = AddMetal( CSX, 'metal' ); %create a PEC called 'metal'
%   start = [0       -width/2  height];
%   stop  = [length  +width/2  0     ];
%   [CSX,port] = AddMSLPort( CSX, 0, 1, 'metal', start, stop, 'x', [0 0 -1], ...
%                           'ExcitePort', true, 'Feed_R', 50 )
% Explanation:
%   - this defines a MSL in x-direction (dir='x')
%     --> the wave travels along the x-direction
%   - with an e-field excitation in -z-direction (evec=[0 0 -1])
%   - the excitation is active and placed at x=start(1) ('ExcitePort', true)
%   - a 50 Ohm lumped port resistance is placed at x=start(1) ('Feed_R', 50)
%   - the width-direction is determined by the cross product of the
%       direction of propagation (dir='x') and the excitation vector
%       (evec=[0 0 -1]), in this case it is the y-direction
%   - the MSL-metal is created in a xy-plane at a height at z=start(3)
%     --> It is important to define the MSL height in the start coordinate!
%   - the ground (xy-plane, not defined by the port) is assumed at z=stop(3)
%     --> The reference plane (ground) is defined in the stop coordinate!
%
% Sebastian Held <sebastian.held@gmx.de> May 13 2010
% Thorsten Liebig <thorsten.liebig@gmx.de> (c) 2011-2013
%
% See also InitCSX DefineRectGrid AddMetal AddMaterial AddExcitation calcPort

%% validate arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%check mesh
if ~isfield(CSX,'RectilinearGrid')
    error 'mesh needs to be defined! Use DefineRectGrid() first!';
end
if (~isfield(CSX.RectilinearGrid,'XLines') || ~isfield(CSX.RectilinearGrid,'YLines') || ~isfield(CSX.RectilinearGrid,'ZLines'))
    error 'mesh needs to be defined! Use DefineRectGrid() first!';
end

% check dir
dir = DirChar2Int(dir);

% check evec
if ~(evec(1) == evec(2) == 0) && ~(evec(1) == evec(3) == 0) && ~(evec(2) == evec(3) == 0) || (sum(evec) == 0)
	error 'evec must have exactly one component ~= 0'
end
evec0 = evec ./ sum(evec); % evec0 is a unit vector

%set defaults
feed_shift = 0;
feed_R = inf; %(default is open, no resistance)
excite = false;
measplanepos = nan;
PortNamePrefix = '';

excite_args = {};

%% read optional arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n=1:2:numel(varargin)
    if (strcmp(varargin{n},'FeedShift')==1);
        feed_shift = varargin{n+1};
        if (numel(feed_shift)>1)
            error 'FeedShift must be a scalar value'
        end
    elseif (strcmp(varargin{n},'Feed_R')==1);
        feed_R = varargin{n+1};
        if (numel(feed_R)>1)
            error 'Feed_R must be a scalar value'
        end
    elseif (strcmp(varargin{n},'MeasPlaneShift')==1);
        measplanepos = varargin{n+1};
        if (numel(measplanepos)>1)
            error 'MeasPlaneShift must be a scalar value'
        end
    elseif (strcmp(varargin{n},'ExcitePort')==1);
        if ischar(varargin{n+1})
            warning('CSXCAD:AddMSLPort','deprecated: a string as excite option is no longer supported and will be removed in the future, please use true or false');
            if ~isempty(excite)
                excite = true;
            else
                excite = false;
            end
        else
            excite = varargin{n+1};
        end
    elseif (strcmpi(varargin{n},'PortNamePrefix'))
        PortNamePrefix = varargin{n+1};
    else
        excite_args{end+1} = varargin{n};
        excite_args{end+1} = varargin{n+1};
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% normalize start and stop
nstart = min( [start;stop] );
nstop  = max( [start;stop] );

% determine index (1, 2 or 3) of propagation (length of MSL)
idx_prop = dir + 1;

% determine index (1, 2 or 3) of width of MSL
dir = [0 0 0];
dir(idx_prop) = 1;
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
CSX = AddBox( CSX, materialname, prio, MSL_start, MSL_stop );

if isnan(measplanepos)
    measplanepos = (nstart(idx_prop)+nstop(idx_prop))/2;
else
    measplanepos = start(idx_prop)+direction*measplanepos;
end

% calculate position of the voltage probes
try
	mesh{1} = sort(CSX.RectilinearGrid.XLines);
	mesh{2} = sort(CSX.RectilinearGrid.YLines);
	mesh{3} = sort(CSX.RectilinearGrid.ZLines);
	meshlines = interp1( mesh{idx_prop}, 1:numel(mesh{idx_prop}), measplanepos, 'nearest' );
	meshlines = mesh{idx_prop}(meshlines-1:meshlines+1); % get three lines (approx. at center)
	if direction == -1
		meshlines = fliplr(meshlines);
	end
	MSL_w2 = interp1( mesh{idx_width}, 1:numel(mesh{idx_width}), (nstart(idx_width)+nstop(idx_width))/2, 'nearest' );
	MSL_w2 = mesh{idx_width}(MSL_w2); % get e-line at center of MSL (MSL_width/2)
	v1_start(idx_prop)   = meshlines(1);
	v1_start(idx_width)  = MSL_w2;
	v1_start(idx_height) = start(idx_height);
	v1_stop  = v1_start;
	v1_stop(idx_height)  = stop(idx_height);
	v2_start = v1_start;
	v2_stop  = v1_stop;
	v2_start(idx_prop)   = meshlines(2);
	v2_stop(idx_prop)    = meshlines(2);
	v3_start = v2_start;
	v3_stop  = v2_stop;
	v3_start(idx_prop)   = meshlines(3);
	v3_stop(idx_prop)    = meshlines(3);
catch
	error('Unable to place voltage probe on mesh; check the location of the MSL and the probe (MeasPlaneShift), and make sure that the mesh is large enough');
end

% calculate position of the current probes
try
	idx = interp1( mesh{idx_width}, 1:numel(mesh{idx_width}), nstart(idx_width), 'nearest' );
	i1_start(idx_width)  = mesh{idx_width}(idx) - diff(mesh{idx_width}(idx-1:idx))/2;
	idx = interp1( mesh{idx_height}, 1:numel(mesh{idx_height}), start(idx_height), 'nearest' );
	i1_start(idx_height) = mesh{idx_height}(idx-1) - diff(mesh{idx_height}(idx-2:idx-1))/2;
	i1_stop(idx_height)  = mesh{idx_height}(idx+1) + diff(mesh{idx_height}(idx+1:idx+2))/2;
	i1_start(idx_prop)   = sum(meshlines(1:2))/2;
	i1_stop(idx_prop)    = i1_start(idx_prop);
	idx = interp1( mesh{idx_width}, 1:numel(mesh{idx_width}), nstop(idx_width), 'nearest' );
	i1_stop(idx_width)   = mesh{idx_width}(idx) + diff(mesh{idx_width}(idx:idx+1))/2;
	i2_start = i1_start;
	i2_stop  = i1_stop;
	i2_start(idx_prop)   = sum(meshlines(2:3))/2;
	i2_stop(idx_prop)    = i2_start(idx_prop);
catch
  error('Unable to place current probe on mesh; check the location of the MSL, and make sure that the mesh is large enough');
end

% create the probes
port.U_filename{1} = [PortNamePrefix 'port_ut' num2str(portnr) 'A'];
% weight = sign(stop(idx_height)-start(idx_height))
weight = 1;
CSX = AddProbe( CSX, port.U_filename{1}, 0, 'weight', weight );
CSX = AddBox( CSX, port.U_filename{1}, prio, v1_start, v1_stop );
port.U_filename{2} = [PortNamePrefix 'port_ut' num2str(portnr) 'B'];
CSX = AddProbe( CSX, port.U_filename{2}, 0, 'weight', weight );
CSX = AddBox( CSX, port.U_filename{2}, prio, v2_start, v2_stop );
port.U_filename{3} = [PortNamePrefix 'port_ut' num2str(portnr) 'C'];
CSX = AddProbe( CSX, port.U_filename{3}, 0, 'weight', weight );
CSX = AddBox( CSX, port.U_filename{3}, prio, v3_start, v3_stop );

weight = direction;
port.I_filename{1} = [PortNamePrefix 'port_it' num2str(portnr) 'A'];
CSX = AddProbe( CSX, port.I_filename{1}, 1, 'weight', weight );
CSX = AddBox( CSX, port.I_filename{1}, prio, i1_start, i1_stop );
port.I_filename{2} = [PortNamePrefix 'port_it' num2str(portnr) 'B'];
CSX = AddProbe( CSX, port.I_filename{2}, 1,'weight', weight );
CSX = AddBox( CSX, port.I_filename{2}, prio, i2_start, i2_stop );

% create port structure
port.LengthScale = 1;
if ((CSX.ATTRIBUTE.CoordSystem==1) && (idx_prop==2))
    port.LengthScale = MSL_stop(idx_height);
end
port.nr = portnr;
port.type = 'MSL';
port.drawingunit = CSX.RectilinearGrid.ATTRIBUTE.DeltaUnit;
port.v_delta = diff(meshlines)*port.LengthScale;
port.i_delta = diff( meshlines(1:end-1) + diff(meshlines)/2 )*port.LengthScale;
port.direction = direction;
port.excite = 0;
port.measplanepos = abs(v2_start(idx_prop) - start(idx_prop))*port.LengthScale;
% port

% create excitation (if enabled) and port resistance
try
	meshline = interp1( mesh{idx_prop}, 1:numel(mesh{idx_prop}), start(idx_prop) + feed_shift*direction, 'nearest' );
	ex_start(idx_prop)   = mesh{idx_prop}(meshline) ;
	ex_start(idx_width)  = nstart(idx_width);
	ex_start(idx_height) = nstart(idx_height);
	ex_stop(idx_prop)    = ex_start(idx_prop);
	ex_stop(idx_width)   = nstop(idx_width);
	ex_stop(idx_height)  = nstop(idx_height);
catch
  error('Unable to place excitation on mesh; check the location of the MSL and the excitation (FeedShift), and make sure that the mesh is large enough');
end

port.excite = 0;
if excite
    port.excite = 1;
    CSX = AddExcitation( CSX, [PortNamePrefix 'port_excite_' num2str(portnr)], 0, evec, excite_args{:} );
    CSX = AddBox( CSX, [PortNamePrefix 'port_excite_' num2str(portnr)], prio, ex_start, ex_stop );
end

%% MSL resistance at start of MSL line
ex_start(idx_prop) = start(idx_prop);
ex_stop(idx_prop) = ex_start(idx_prop);

if (feed_R > 0) && ~isinf(feed_R)
    CSX = AddLumpedElement( CSX, [PortNamePrefix 'port_resist_' int2str(portnr)], idx_height-1, 'R', feed_R );
    CSX = AddBox( CSX, [PortNamePrefix 'port_resist_' int2str(portnr)], prio, ex_start, ex_stop );
elseif isinf(feed_R)
    % do nothing --> open port
elseif feed_R == 0
    %port "resistance" as metal
    CSX = AddBox( CSX, materialname, prio, ex_start, ex_stop );
else
    error('openEMS:AddMSLPort','MSL port with resistance <= 0 it not possible');
end
end
