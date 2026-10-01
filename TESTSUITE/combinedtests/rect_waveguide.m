function pass = rect_waveguide(varargin)
%pass = rect_waveguide(<key>, <value>, ...)
%
% WR-42 rectangular waveguide with a TE10 mode matching port at each end. The
% ports excite and measure through the analytic mode profile (probe types
% 10/11), so a port that launches a clean mode into a well matched PML shows
% almost no reflection and no loss.
%
% The phase of S21 is checked against the analytic waveguide dispersion
%
%     beta = sqrt(k0^2 - kc^2),   kc = pi/a  for TE10
%
% which is the part of this that the solver has to get right -- a simple
% beta = k0 would be off by a third at 24 GHz in a WR-42.
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
physical_constants;

% LIMITS
limit_S11_dB  = -20;   % matched mode matching port into a PML
limit_S21_dB  = 0.5;   % lossless waveguide
limit_phase   = 0.02;  % 2% on the S21 phase delay

% geometry, in mm: WR-42
unit     = 1e-3;
a        = 10.7;
b        =  4.3;
wg_length = 50;

f_start = 20e9;
f_stop  = 26e9;
f       = linspace(f_start, f_stop, 201);

kc = pi/(a*unit);        % TE10 cutoff wavenumber
fc = c0*kc/(2*pi);

Sim_Path = ts_sim_path(mfilename('fullpath'));
Sim_CSX  = 'rect_wg.xml';

%% FDTD setup
FDTD = InitFDTD('NrTS', 20000, 'EndCriteria', 1e-5);
FDTD = SetGaussExcite(FDTD, 0.5*(f_start+f_stop), 0.5*(f_stop-f_start));
FDTD = SetBoundaryCond(FDTD, {'PEC','PEC','PEC','PEC','PML_8','PML_8'});

%% mesh
% uniform, about lambda/30 at the top of the band: a waveguide needs no
% grading, and the port faces then sit exactly on mesh lines
CSX = InitCSX();
dz = wg_length / ceil(wg_length / (c0/f_stop/unit/30));
mesh.x = linspace(0, a, ceil(a/dz)+1);
mesh.y = linspace(0, b, ceil(b/dz)+1);
mesh.z = 0:dz:wg_length;
CSX = DefineRectGrid(CSX, unit, mesh);

% the ports sit ten cells away from the absorbing boundaries
z_p1 = [10 15] * dz;
z_p2 = wg_length - [10 15] * dz;

%% TE10 mode matching ports
[CSX,port{1}] = AddRectWaveGuidePort(CSX, 0, 1, [0 0 z_p1(1)], [a b z_p1(2)], ...
                                     'z', a*unit, b*unit, 'TE10', 1);
[CSX,port{2}] = AddRectWaveGuidePort(CSX, 0, 2, [0 0 z_p2(1)], [a b z_p2(2)], ...
                                     'z', a*unit, b*unit, 'TE10', 0);

%% run
WriteOpenEMS([Sim_Path '/' Sim_CSX], FDTD, CSX);
Settings.LogFile = [Sim_Path '/openEMS.log'];
Settings.Silent  = opt.Silent;
RunOpenEMS(Sim_Path, Sim_CSX, opt.openEMS_opts, Settings);

%% analysis
port = calcPort(port, Sim_Path, f);

S11 = port{1}.uf.ref ./ port{1}.uf.inc;
S21 = port{2}.uf.ref ./ port{1}.uf.inc;

k    = 2*pi*f/c0;
beta = sqrt(k.^2 - kc^2);
% distance between the two measurement planes, which snap to mesh lines
meas_distance = abs(port{2}.measplanepos - port{1}.measplanepos) * unit;

if opt.Plots
    figure
    plot(f/1e9, 20*log10(abs([S11; S21])), 'LineWidth', 2);
    xlabel('Frequency (GHz)'); ylabel('S-parameter (dB)');
    legend({'S11','S21'});
    title(sprintf('WR-42, TE10 cutoff at %.2f GHz', fc/1e9));
end

pass = ts_check('TE10 port matched', ...
                max(20*log10(abs(S11))) <= limit_S11_dB, ...
                'max. S11 = %.1f dB (limit %.1f dB)', ...
                max(20*log10(abs(S11))), limit_S11_dB);

pass = ts_check('lossless transmission', ...
                max(abs(20*log10(abs(S21)))) <= limit_S21_dB, ...
                'S21 = %.2f .. %.2f dB (limit +/-%.2f dB)', ...
                min(20*log10(abs(S21))), max(20*log10(abs(S21))), limit_S21_dB) && pass;

% total phase delay between the measurement planes; the measured phase is
% lifted onto the same branch, the dispersion is what is compared
phase_ref  = -beta * meas_distance;
phase_meas = unwrap(angle(S21));
phase_meas = phase_meas + 2*pi*round((phase_ref(1) - phase_meas(1))/(2*pi));
pass = ts_check_rel('S21 phase follows the waveguide dispersion', ...
                    phase_meas, phase_ref, limit_phase) && pass;

pass = ts_finish(opt, mfilename, pass, Sim_Path);
