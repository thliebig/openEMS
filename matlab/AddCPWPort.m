function [CSX,port] = AddCPWPort( CSX, prio, portnr, materialname, start, stop, gap_width, dir, evec, varargin )
% [CSX,port] = AddCPWPort( CSX, prio, portnr, materialname, start, stop, gap_width, dir, evec, varargin )
%
% CSX:          CSX-object created by InitCSX()
% prio:         priority for excitation and probe boxes
% portnr:       (integer) number of the port
% materialname: property for the CPW (created by AddMetal())
% start:        3D start rowvector for port definition
% stop:         3D end rowvector for port definition
% gap_width:    width of the CPW gap (left and right)
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
%   start = [0       -width/2  0];
%   stop  = [length  +width/2  0];
%   [CSX,port] = AddCPWPort( CSX, 0, 1, 'metal', start, stop, gap_width, 'x', ...
%                                   [0 0 -1], 'ExcitePort', true, 'Feed_R', 50 )
% Explanation:
%   - this defines a stripline in x-direction (dir='x')
%     --> the wave travels along the x-direction
%   - with an e-field excitation in -z-direction (evec=[0 0 -1])
%   - the excitation is active and placed at x=start(1) ('ExcitePort', true)
%   - a 50 Ohm lumped port resistance is placed at x=start(1) ('Feed_R', 50)
%   - the width-direction is determined by the cross product of the
%       direction of propagtion (dir='x') and the excitation vector
%       (evec=[0 0 -1]), in this case it is the y-direction
%   - the stripline-metal is created in a xy-plane at a height at z=start(3)
%     --> The upper and lower reference plane (ground) must be defined by
%     the user
%
% Thorsten Liebig <thorsten.liebig@gmx.de> (c) 2014
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
            warning('CSXCAD:AddCPWPort','depreceated: a string as excite option is no longer supported and will be removed in the future, please use true or false');
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

% determine index (1, 2 or 3) of propagation (length of CPW)
idx_prop = dir + 1;

% determine index (1, 2 or 3) of width of CPW
dir = [0 0 0];
dir(idx_prop) = 1;
idx_height = abs(cross(dir,evec0)) * [1;2;3];

% determine index (1, 2 or 3) of height
idx_width = abs(evec0) * [1;2;3];


if (start(idx_height)~=stop(idx_height))
    error('openEMS:AddCPWPort','start/stop in height direction must be equal');
end

% direction of propagation
if stop(idx_prop)-start(idx_prop) > 0
    direction = +1;
else
    direction = -1;
end

% create the metal/material for the CPW
SL_start = start;
SL_stop = stop;
CSX = AddBox( CSX, materialname, prio, SL_start, SL_stop );

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
	SL_w2 = interp1( mesh{idx_width}, 1:numel(mesh{idx_width}), (nstart(idx_width)+nstop(idx_width))/2, 'nearest' );
	SL_w2 = mesh{idx_width}(SL_w2); % get e-line at center of CPW (SL_width/2)
	v1_start(idx_prop)   = meshlines(1);
	v1_start(idx_width)  = (nstart(idx_width)+nstop(idx_width))/2;
	v1_start(idx_height) = start(idx_height);
	v1_stop  = v1_start;
	v2_start = v1_start;
	v2_stop  = v1_stop;
	v2_start(idx_prop)   = meshlines(2);
	v2_stop(idx_prop)    = meshlines(2);
	v3_start = v2_start;
	v3_stop  = v2_stop;
	v3_start(idx_prop)   = meshlines(3);
	v3_stop(idx_prop)    = meshlines(3);
catch
	error('Unable to place voltage probe on mesh; check the location of the CPW and the probe (MeasPlaneShift), and make sure that the mesh is large enough');
end

width_add_start = [0 0 0];
width_add_stop  = [0 0 0];
width_add_start(idx_width) = (nstop(idx_width)-nstart(idx_width))/2;
width_add_stop(idx_width)  = (nstop(idx_width)-nstart(idx_width))/2+gap_width;

weight = 0.5;
% create the voltage-probes
port.U_filename{1,1} = [PortNamePrefix 'port_ut' num2str(portnr) 'A1'];
CSX = AddProbe( CSX, port.U_filename{1,1}, 0, 'weight', weight );
CSX = AddBox( CSX, port.U_filename{1,1}, prio, v1_start-width_add_start, v1_stop-width_add_stop);

port.U_filename{1,2} = [PortNamePrefix 'port_ut' num2str(portnr) 'A2'];
CSX = AddProbe( CSX, port.U_filename{1,2}, 0, 'weight', weight );
CSX = AddBox( CSX, port.U_filename{1,2}, prio, v1_start+width_add_start, v1_stop+width_add_stop);


port.U_filename{2,1} = [PortNamePrefix 'port_ut' num2str(portnr) 'B1'];
CSX = AddProbe( CSX, port.U_filename{2,1}, 0, 'weight', weight );
CSX = AddBox( CSX, port.U_filename{2,1}, prio, v2_start-width_add_start, v2_stop-width_add_stop );

port.U_filename{2,2} = [PortNamePrefix 'port_ut' num2str(portnr) 'B2'];
CSX = AddProbe( CSX, port.U_filename{2,2}, 0, 'weight', weight );
CSX = AddBox( CSX, port.U_filename{2,2}, prio, v2_start+width_add_start, v2_stop+width_add_stop );


port.U_filename{3,1} = [PortNamePrefix 'port_ut' num2str(portnr) 'C1'];
CSX = AddProbe( CSX, port.U_filename{3,1}, 0, 'weight', weight );
CSX = AddBox( CSX, port.U_filename{3,1}, prio, v3_start-width_add_start, v3_stop-width_add_stop );

port.U_filename{3,2} = [PortNamePrefix 'port_ut' num2str(portnr) 'C2'];
CSX = AddProbe( CSX, port.U_filename{3,2}, 0, 'weight', weight );
CSX = AddBox( CSX, port.U_filename{3,2}, prio, v3_start+width_add_start, v3_stop+width_add_stop );

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
	error('Unable to place current probe on mesh; check the location of the CPW and the probe (MeasPlaneShift), and make sure that the mesh is large enough');
end

% create the curr-probes
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
    port.LengthScale = SL_stop(idx_height);
end
port.nr = portnr;
port.type = 'CPW';
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
	ex_start(idx_width)  = (nstart(idx_width)+nstop(idx_width))/2;
	ex_start(idx_height) = nstart(idx_height);
	ex_stop(idx_prop)    = ex_start(idx_prop);
	ex_stop(idx_width)   = (nstart(idx_width)+nstop(idx_width))/2;
	ex_stop(idx_height)  = nstop(idx_height);
catch
	error('Unable to place excitation on mesh; check the location of the CPW and the excitation (FeedShift), and make sure that the mesh is large enough');
end

port.excite = 0;
if excite
    port.excite = 1;
    CSX = AddExcitation( CSX, [PortNamePrefix 'port_excite_1_' num2str(portnr)], 0, evec, excite_args{:} );
    CSX = AddBox( CSX, [PortNamePrefix 'port_excite_1_' num2str(portnr)], prio, ex_start-width_add_start, ex_stop-width_add_stop );
    CSX = AddExcitation( CSX, [PortNamePrefix 'port_excite_2_' num2str(portnr)], 0, -evec, excite_args{:} );
    CSX = AddBox( CSX, [PortNamePrefix 'port_excite_2_' num2str(portnr)], prio, ex_start+width_add_start, ex_stop+width_add_stop );
end

%% CPW resistance at start of CPW line
ex_start(idx_prop) = start(idx_prop);
ex_stop(idx_prop) = ex_start(idx_prop);

if (feed_R > 0) && ~isinf(feed_R)
    CSX = AddLumpedElement( CSX, [PortNamePrefix 'port_resist_' int2str(portnr)], idx_width-1, 'R', 2*feed_R );
    CSX = AddBox( CSX, [PortNamePrefix 'port_resist_' int2str(portnr)], prio, ex_start-width_add_start, ex_stop-width_add_stop );
    CSX = AddBox( CSX, [PortNamePrefix 'port_resist_' int2str(portnr)], prio, ex_start+width_add_start, ex_stop+width_add_stop );
elseif isinf(feed_R)
    % do nothing --> open port
elseif feed_R == 0
    %port "resistance" as metal
    CSX = AddBox( CSX, materialname, prio, ex_start-width_add_start, ex_stop-width_add_stop );
    CSX = AddBox( CSX, materialname, prio, ex_start+width_add_start, ex_stop+width_add_stop );
else
    error('openEMS:AddCPWPort','CPW port with resistance <= 0 it not possible');
end
end
