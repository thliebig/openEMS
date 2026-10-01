function pass = coax_cylindrical(varargin)
%pass = coax_cylindrical(<key>, <value>, ...)
%
% Coaxial line in cylindrical coordinates, excited with the analytic 1/rho TEM
% profile. This is what covers Operator_Cylinder and the closed 2*pi alpha
% mesh; the cartesian Coax test cannot.
%
% Checked against closed form results of the ideal line:
%
%   - the mode is TEM, so beta = k0
%   - E_rho ~ 1/rho, so the voltage between r_i and any radius r follows
%     log(r/r_i), independent of the mesh
%   - the TEM mode is rotationally symmetric, so the voltage does not depend
%     on alpha at all -- that is what breaks first if the duplicated mesh line
%     of a closed alpha range is mishandled
%
% Not checked here: the characteristic impedance. It needs a current probe
% around the inner conductor, and such a probe reports one alpha cell too
% little on a closed cylindrical mesh (ZL comes out a factor N/(N-1) high for
% N alpha cells, measured over N = 8..64). Add the check once that is fixed.
%
% See ts_options for the accepted options.
%
% openEMS testsuite
% -----------------
%
% See also run_testsuite, ts_options, Coax

% the shared check helpers, also when this test is run on its own
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'helperscripts'));

opt = ts_options(varargin{:});
physical_constants;

% LIMITS
limit_beta     = 0.01;   % 1% on the propagation constant
limit_profile  = 0.01;   % 1% on the radial voltage distribution
limit_symmetry = 1e-9;   % the rotational symmetry is exact, not approximate

% geometry, in mm
unit        = 1e-3;
coax_rad_i  = 100;
coax_rad_o  = 230;
coax_rad_p  = 160;      % a mesh line in between, for the radial profile
coax_length = 1000;
N_alpha     = 17;       % the TEM mode is rotationally symmetric, few lines do

f0 = 0.5e9;
f  = linspace(0.1e9, 1e9, 181);

Sim_Path = ts_sim_path(mfilename('fullpath'));
Sim_CSX  = 'coax_cyl.xml';

%% FDTD setup
FDTD = InitFDTD('CoordSystem', 1, 'NrTS', 200000, 'EndCriteria', 1e-5);
FDTD = SetGaussExcite(FDTD, f0, f0);
FDTD = SetBoundaryCond(FDTD, {'PEC','PEC','PEC','PEC','PEC','PML_8'});

%% mesh: rho, alpha, z
CSX = InitCSX('CoordSystem', 1);
mesh.x = coax_rad_i : 10 : coax_rad_o;
mesh.y = linspace(0, 2*pi, N_alpha);
mesh.z = 0 : 10 : coax_length;
CSX = DefineRectGrid(CSX, unit, mesh);

%% excitation: the TEM profile E_rho ~ 1/rho over the full cross-section
CSX = AddExcitation(CSX, 'excite', 0, [1 0 0]);
CSX = SetExcitationWeight(CSX, 'excite', {'1/rho', 0, 0});
CSX = AddBox(CSX, 'excite', 0, [mesh.x(1) mesh.y(1) 0], [mesh.x(end) mesh.y(end) 0]);

%% voltage probes: radial lines
z1 = mesh.z(11);
z2 = mesh.z(61);
CSX = AddProbe(CSX, 'u1', 0);           % full gap, alpha = 0
CSX = AddBox(CSX, 'u1', 0, [coax_rad_i 0 z1], [coax_rad_o 0 z1]);
CSX = AddProbe(CSX, 'u1_rot', 0);       % full gap, a quarter turn away
CSX = AddBox(CSX, 'u1_rot', 0, [coax_rad_i pi/2 z1], [coax_rad_o pi/2 z1]);
CSX = AddProbe(CSX, 'u1_part', 0);      % inner part of the gap only
CSX = AddBox(CSX, 'u1_part', 0, [coax_rad_i 0 z1], [coax_rad_p 0 z1]);
CSX = AddProbe(CSX, 'u2', 0);           % full gap, further down the line
CSX = AddBox(CSX, 'u2', 0, [coax_rad_i 0 z2], [coax_rad_o 0 z2]);

%% run
WriteOpenEMS([Sim_Path '/' Sim_CSX], FDTD, CSX);
Settings.LogFile = [Sim_Path '/openEMS.log'];
Settings.Silent  = opt.Silent;
RunOpenEMS(Sim_Path, Sim_CSX, opt.openEMS_opts, Settings);

%% analysis
U = ReadUI({'u1','u1_rot','u1_part','u2'}, Sim_Path, f);
u1      = U.FD{1}.val;
u1_rot  = U.FD{2}.val;
u1_part = U.FD{3}.val;
u2      = U.FD{4}.val;

% the line ends in a PML, so u2/u1 is the phase delay of a forward travelling
% wave over the distance between the two probes
d          = (z2 - z1) * unit;
k0         = 2*pi*f/c0;
phase_ref  = -k0 * d;
phase_meas = unwrap(angle(u2 ./ u1));
phase_meas = phase_meas + 2*pi*round((phase_ref(1) - phase_meas(1))/(2*pi));

profile     = abs(u1_part ./ u1);
profile_ref = log(coax_rad_p/coax_rad_i) / log(coax_rad_o/coax_rad_i);

if opt.Plots
    figure
    plot(f/1e9, [phase_meas; phase_ref], 'LineWidth', 2);
    xlabel('Frequency (GHz)'); ylabel('phase of u_2/u_1 (rad)');
    legend({'simulated','-k_0 d'});

    figure
    plot(f/1e9, profile, 'k-', 'LineWidth', 2);
    hold on
    plot(f/1e9, profile_ref*ones(size(f)), 'r--', 'LineWidth', 2);
    xlabel('Frequency (GHz)'); ylabel('u(r_i..r_p) / u(r_i..r_o)');
    legend({'simulated','log(r_p/r_i)/log(r_o/r_i)'});
end

pass = ts_check_rel('TEM propagation, beta = k0', phase_meas, phase_ref, limit_beta);

pass = ts_check_rel('radial voltage distribution follows log(r)', ...
                    profile, profile_ref, limit_profile) && pass;

sym_err = max(abs(u1_rot./u1 - 1));
pass = ts_check('TEM mode is rotationally symmetric', sym_err <= limit_symmetry, ...
                'max. deviation between alpha = 0 and 90 deg: %.3g (limit %.3g)', ...
                sym_err, limit_symmetry) && pass;

pass = ts_finish(opt, mfilename, pass, Sim_Path);
