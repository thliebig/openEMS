function pass = port_setup(varargin)
%pass = port_setup(<key>, <value>, ...)
%
% Check what the port helpers put into the CSX structure, without running a
% simulation. What is checked here is the measurement definition itself: the
% probe geometry, the probe weights and the sign conventions that every
% S-parameter in this suite is built on, plus the cutoff wavenumbers of the
% rectangular waveguide modes.
%
% A sign or a weight that silently changes here shows up as a wrong
% S-parameter in a full simulation, where it is much harder to pin down.
%
% See ts_options for the accepted options.
%
% openEMS testsuite
% -----------------
%
% See also run_testsuite, ts_options

% the shared check helpers, also when this test is run on its own
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'helperscripts'));

opt = ts_options(varargin{:});
pass = 1;

unit = 1e-3;

%% lumped port %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a 50 Ohm port across a 2 mm gap, excited in +z
start = [-1 -2  0];
stop  = [ 1  2 10];
CSX = base_csx(unit);
[CSX, port] = AddLumpedPort(CSX, 10, 1, 50, start, stop, [0 0 1], true);

pass = ts_check('lumped port direction', port.direction == +1, ...
                'stop above start in z gives +1, got %+d', port.direction) && pass;

le = find_prop(CSX, 'LumpedElement', 'port_resist_1');
pass = ts_check('lumped port resistor', ...
                ~isempty(le) && le.ATTRIBUTE.R == 50 && le.ATTRIBUTE.Direction == 2, ...
                'R = %g Ohm in direction %d (z)', le.ATTRIBUTE.R, le.ATTRIBUTE.Direction) && pass;

% the voltage probe integrates along the port direction through the centre of
% the gap, the current probe encloses the port in the plane normal to it
up = find_prop(CSX, 'ProbeBox', port.U_filename);
ip = find_prop(CSX, 'ProbeBox', port.I_filename);
u_box = [up.Primitives.Box{1}.P1.ATTRIBUTE; up.Primitives.Box{1}.P2.ATTRIBUTE];
i_box = [ip.Primitives.Box{1}.P1.ATTRIBUTE; ip.Primitives.Box{1}.P2.ATTRIBUTE];
u_1 = [u_box(1).X u_box(1).Y u_box(1).Z];
u_2 = [u_box(2).X u_box(2).Y u_box(2).Z];
i_1 = [i_box(1).X i_box(1).Y i_box(1).Z];
i_2 = [i_box(2).X i_box(2).Y i_box(2).Z];

pass = ts_check('lumped port voltage probe geometry', ...
                isequal(u_1, [0 0 0]) && isequal(u_2, [0 0 10]), ...
                'line (%g,%g,%g)..(%g,%g,%g), expected the gap centre line', ...
                u_1, u_2) && pass;
pass = ts_check('lumped port current probe geometry', ...
                isequal(i_1, [-1 -2 5]) && isequal(i_2, [1 2 5]), ...
                'plane (%g,%g,%g)..(%g,%g,%g), expected the mid-gap cross-section', ...
                i_1, i_2) && pass;

% u is measured against the port direction, i along it: that is the
% convention calcLumpedPort assumes for the incident/reflected split
pass = ts_check('lumped port probe weights', ...
                up.ATTRIBUTE.Weight == -1 && ip.ATTRIBUTE.Weight == +1 && ...
                ip.ATTRIBUTE.NormDir == 2, ...
                'U weight %+g, I weight %+g, I normal direction %d', ...
                up.ATTRIBUTE.Weight, ip.ATTRIBUTE.Weight, ip.ATTRIBUTE.NormDir) && pass;

pass = ts_check('lumped port excitation', ...
                ~isempty(find_prop(CSX, 'Excitation', 'port_excite_1')), ...
                'excite=true adds the excitation property') && pass;

% a port pointing the other way flips the current probe weight
CSX = base_csx(unit);
[CSX, port_rev] = AddLumpedPort(CSX, 10, 1, 50, stop, start, [0 0 1], false);
ip_rev = find_prop(CSX, 'ProbeBox', port_rev.I_filename);
pass = ts_check('reversed lumped port', ...
                port_rev.direction == -1 && ip_rev.ATTRIBUTE.Weight == -1 && ...
                isempty(find_prop(CSX, 'Excitation', 'port_excite_1')), ...
                'direction %+d, I weight %+g, no excitation for excite=false', ...
                port_rev.direction, ip_rev.ATTRIBUTE.Weight) && pass;

%% rectangular waveguide port %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WR-42, cutoff wavenumbers after Pozar, Microwave Engineering
a = 10.7e-3;
b =  4.3e-3;
modes = {'TE10','TE11','TE20','TE21'};
mn    = [1 0; 1 1; 2 0; 2 1];

kc     = zeros(1,numel(modes));
kc_ref = zeros(1,numel(modes));
for n=1:numel(modes)
    CSX = base_csx(unit);
    [CSX, wg] = AddRectWaveGuidePort(CSX, 0, 1, [0 0 0], [a/unit b/unit 20], ...
                                     'z', a, b, modes{n}, 0);
    kc(n)     = wg.kc;
    kc_ref(n) = sqrt((mn(n,1)*pi/a)^2 + (mn(n,2)*pi/b)^2);
end
pass = ts_check_rel('rect. waveguide cutoff wavenumbers', kc, kc_ref, 1e-12) && pass;

% the mode profile has to be written in the two coordinates transverse to the
% direction of propagation, and the mode probes have to be mode matching ones
for dir = {'x','y','z'}
    CSX = base_csx(unit);
    [CSX, wg] = AddRectWaveGuidePort(CSX, 0, 1, [0 0 0], [20 20 20], ...
                                     dir{1}, a, b, 'TE10', 0);
    mode = find_prop(CSX, 'ProbeBox', wg.U_filename).ModeFunction.ATTRIBUTE;
    % the component along the direction of propagation is written as 0
    func = strtrim([as_str(mode.X) ' ' as_str(mode.Y) ' ' as_str(mode.Z)]);
    pass = ts_check(['rect. waveguide TE10 mode profile, ' dir{1} '-direction'], ...
                    isempty(strfind(func, dir{1})), ...
                    'transverse coordinates only, got "%s"', func) && pass;
end

CSX = base_csx(unit);
[CSX, wg] = AddRectWaveGuidePort(CSX, 0, 1, [0 0 0], [20 20 20], 'z', a, b, 'TE10', 1);
up = find_prop(CSX, 'ProbeBox', wg.U_filename);
ip = find_prop(CSX, 'ProbeBox', wg.I_filename);
pass = ts_check('rect. waveguide mode matching probes', ...
                up.ATTRIBUTE.Type == 10 && ip.ATTRIBUTE.Type == 11, ...
                'probe types %d and %d (10/11 = waveguide mode matching)', ...
                up.ATTRIBUTE.Type, ip.ATTRIBUTE.Type) && pass;

pass = ts_finish(opt, mfilename, pass);

end


function CSX = base_csx(unit)
% a mesh big enough for every port built here; the port helpers need one
CSX = InitCSX();
mesh.x = linspace(-20, 20, 41);
mesh.y = linspace(-20, 20, 41);
mesh.z = linspace(  0, 40, 41);
CSX = DefineRectGrid(CSX, unit, mesh);
end


function prop = find_prop(CSX, type, name)
% the CSX property of the given type and name, [] if there is none
prop = [];
if ~isfield(CSX,'Properties') || ~isfield(CSX.Properties, type)
    return
end
for n=1:numel(CSX.Properties.(type))
    if strcmp(CSX.Properties.(type){n}.ATTRIBUTE.Name, name)
        prop = CSX.Properties.(type){n};
        return
    end
end
end


function s = as_str(v)
% mode function components are either an expression or the number 0
if ischar(v)
    s = v;
else
    s = num2str(v);
end
end
