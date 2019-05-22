function [CSX,port] = AddCoaxialPort( CSX, prio, portnr, pec_name, materialname, start, stop, dir, r_i, r_o, r_os, varargin )
% function [CSX,port] = AddCoaxialPort( CSX, prio, portnr, pec_name, materialname, start, stop, dir, r_i, r_o, r_os, varargin )
%
% CSX:          CSX-object created by InitCSX()
% prio:         priority for excitation and probe boxes
% portnr:       (integer) number of the port
% pec_name:     metal property for coaxial inner/outer conductor (created by AddMetal())
% materialname: substrate property for coaxial line (created by AddMaterial())
%               Note: this may be empty for an "air filled" coaxial line
% start:        3D start rowvector for coaxial cable axis
% stop:         3D end rowvector for coaxial cable axis
% dir:          direction of wave propagation (choices: 0, 1, 2 or 'x','y','z')
% r_i:          inner coaxial radius (in drawing unit)
% r_o:          outer coaxial radius (in drawing unit)
% r_os:         outer shell coaxial radius (in drawing unit)
%
% variable input:
%  varargin:    optional additional excitations options, see also AddExcitation
% 'ExciteAmp'   excitation amplitude of transversal electric field profile,
%               set to 0 (default) for a passive port
% 'FeedShift'   shift to port from start by a given distance in drawing
%               units. Default is 0. Only active if 'ExciteAmp' is set!
% 'Feed_R'      Specify a lumped port resistance. Default is no lumped
%               port resistance --> port has to end in an ABC.
% 'MeasPlaneShift'  Shift the measurement plane from start t a given distance
%               in drawing units. Default is the middle of start/stop.
% 'PortNamePrefix' a prefix to the port name
%
% the mesh must be already initialized
%
% example:
%
% openEMS matlab interface
% -----------------------
% Thorsten Liebig <thorsten.liebig@gmx.de> (c) 2013
%
% See also InitCSX AddMetal AddMaterial AddExcitation calcPort

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

%set defaults
feed_shift = 0;
feed_R = inf; %(default is open, no resitance)
excite_amp = 0;
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
    elseif (strcmp(varargin{n},'ExciteAmp')==1);
        excite_amp = varargin{n+1};
    elseif (strcmpi(varargin{n},'PortNamePrefix'))
        PortNamePrefix = varargin{n+1};
    else
        excite_args{end+1} = varargin{n};
        excite_args{end+1} = varargin{n+1};
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% determine index (1, 2 or 3) of propagation (length of MSL)
idx_prop_n   = dir + 1;
idx_prop_nP  = mod((dir+1),3)+1;
idx_prop_nPP = mod((dir+2),3)+1;

% direction of propagation
if stop(idx_prop_n)-start(idx_prop_n) > 0
    direction = +1;
else
    direction = -1;
end

% create the metal for the coaxial line
CSX = AddCylinder( CSX, pec_name, prio, start, stop, r_i );
CSX = AddCylindricalShell( CSX, pec_name, prio, start, stop, 0.5*(r_o+r_os), r_os-r_o );

% create the material filling for the coaxial line
if (~isempty(materialname))
    CSX = AddCylindricalShell( CSX, materialname, prio-1, start, stop, 0.5*(r_o+r_i), r_o-r_i );
end

if isnan(measplanepos)
    measplanepos = (start(idx_prop_n)+stop(idx_prop_n))/2;
else
    measplanepos = start(idx_prop_n)+direction*measplanepos;
end

% calculate position of the voltage probes
try
	mesh{1} = sort(unique(CSX.RectilinearGrid.XLines));
	mesh{2} = sort(unique(CSX.RectilinearGrid.YLines));
	mesh{3} = sort(unique(CSX.RectilinearGrid.ZLines));
	meshlines = interp1( mesh{idx_prop_n}, 1:numel(mesh{idx_prop_n}), measplanepos, 'nearest' );
	meshlines = mesh{idx_prop_n}(meshlines-1:meshlines+1); % get three lines (approx. at center)
	if direction == -1
		meshlines = fliplr(meshlines);
	end
	v1_start(idx_prop_n)   = meshlines(1);
	v1_start(idx_prop_nP)  = start(idx_prop_nP)+r_i;
	v1_start(idx_prop_nPP) = start(idx_prop_nPP);
	v1_stop  = v1_start;
	v1_stop(idx_prop_nP)  = start(idx_prop_nP)+r_o;
	v2_start = v1_start;
	v2_stop  = v1_stop;
	v2_start(idx_prop_n)   = meshlines(2);
	v2_stop(idx_prop_n)    = meshlines(2);
	v3_start = v2_start;
	v3_stop  = v2_stop;
	v3_start(idx_prop_n)   = meshlines(3);
	v3_stop(idx_prop_n)    = meshlines(3);
catch
	error('Unable to place voltage probe on mesh; check the location of the port and the probe (MeasPlaneShift), and make sure that the mesh is large enough');
end

% calculate position of the current probes
i1_start(idx_prop_n)   = 0.5*(meshlines(1)+meshlines(2));
i1_start(idx_prop_nP)  = start(idx_prop_nP)-r_i-0.1*(r_o-r_i);
i1_start(idx_prop_nPP) = start(idx_prop_nPP)-r_i-0.1*(r_o-r_i);
i1_stop  = i1_start;
i1_stop(idx_prop_nP)  = start(idx_prop_nP)+r_i+0.1*(r_o-r_i);
i1_stop(idx_prop_nPP) = start(idx_prop_nPP)+r_i+0.1*(r_o-r_i);

i2_start = i1_start;
i2_stop  = i1_stop;
i2_start(idx_prop_n)   = 0.5*(meshlines(2)+meshlines(3));
i2_stop(idx_prop_n)    = 0.5*(meshlines(2)+meshlines(3));

% create the probes
port.U_filename{1} = [PortNamePrefix 'port_ut' num2str(portnr) 'A'];
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
port.nr = portnr;
port.type = 'Coaxial';
port.drawingunit = CSX.RectilinearGrid.ATTRIBUTE.DeltaUnit;
port.v_delta = diff(meshlines)*port.LengthScale;
port.i_delta = diff( meshlines(1:end-1) + diff(meshlines)/2 )*port.LengthScale;
port.direction = direction;
port.excite = 0;
port.measplanepos = abs(v2_start(idx_prop_n) - start(idx_prop_n))*port.LengthScale;

port.r_i = r_i;
port.r_o = r_o;

% create excitation (if enabled) and port resistance
try
	meshline = interp1( mesh{idx_prop_n}, 1:numel(mesh{idx_prop_n}), start(idx_prop_n) + feed_shift*direction, 'nearest' );
	min_cell_prop = min(diff(mesh{idx_prop_n}));
	ex_start = start;
	ex_start(idx_prop_n)   = mesh{idx_prop_n}(meshline) - 0.01*min_cell_prop;
	ex_stop  = ex_start;
	ex_stop(idx_prop_n)   = mesh{idx_prop_n}(meshline) + 0.01*min_cell_prop;
catch
	error('Unable to place excitation on mesh; check the location of the port and the excitation (FeedShift), and make sure that the mesh is large enough');
end

port.excite = 0;
if (excite_amp~=0)
    dir_names={'x','y','z'};
    nameX = ['(' dir_names{idx_prop_nP}  '-' num2str(start(idx_prop_nP)) ')'];
    nameY = ['(' dir_names{idx_prop_nPP} '-' num2str(start(idx_prop_nPP)) ')'];

    func_Ex = [ nameX '/(' nameX '*' nameX '+' nameY '*' nameY ') * (sqrt(' nameX '*' nameX '+' nameY '*' nameY ')<' num2str(r_o) ') * (sqrt(' nameX '*' nameX '+' nameY '*' nameY ')>' num2str(r_i) ')'];
    func_Ey = [ nameY '/(' nameX '*' nameX '+' nameY '*' nameY ') * (sqrt(' nameX '*' nameX '+' nameY '*' nameY ')<' num2str(r_o) ') * (sqrt(' nameX '*' nameX '+' nameY '*' nameY ')>' num2str(r_i) ')'];

    func_E{idx_prop_n} = 0;
    func_E{idx_prop_nP} = func_Ex;
    func_E{idx_prop_nPP} = func_Ey;

    port.excite = 1;
    evec = [1 1 1];
    evec(idx_prop_n) = 0;

    CSX = AddExcitation( CSX, [PortNamePrefix 'port_excite_' num2str(portnr)], 0, evec, excite_args{:} );
    CSX = SetExcitationWeight(CSX, [PortNamePrefix 'port_excite_' num2str(portnr)], func_E );
    CSX = AddCylindricalShell(CSX,[PortNamePrefix 'port_excite_' num2str(portnr)],0 ,ex_start,ex_stop,0.5*(r_i+r_o),(r_o-r_i));
end

%% resistance at start of coaxial line
ex_start = start;
ex_stop = stop;
ex_stop(idx_prop_n) = ex_start(idx_prop_n);

if (feed_R > 0) && ~isinf(feed_R)
    error 'feed_R not yet implemented'
elseif isinf(feed_R)
    % do nothing --> open port
elseif feed_R == 0
    %port "resistance" as metal
    CSX = AddBox( CSX, pec_name, prio, ex_start, ex_stop );
    CSX = AddCylindricalShell(CSX, pec_name, prio ,ex_start, ex_stop, 0.5*(r_i+r_o),(r_o-r_i));
else
    error('openEMS:AddCoaxialPort','Coaxial port with resistance <= 0 it not possible');
end
end
