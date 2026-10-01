function pass = msl_port(varargin)
%pass = msl_port(<key>, <value>, ...)
%
% Microstrip line over a PEC ground, air filled, with a transmission line port
% at each end. An air filled microstrip is a true TEM line, so its propagation
% constant is the free space wavenumber and its characteristic impedance is
% frequency independent -- both are analytic references for what calcTLPort
% extracts out of the two voltage and the current probe of a port.
%
% Checked: beta against k0, a real and flat line impedance, and a line that
% neither reflects nor loses power between the two measurement planes.
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
limit_beta       = 0.02;  % 2%, dominated by the numerical dispersion at lambda/30
limit_ZL_imag    = 0.03;  % |imag(ZL)| / real(ZL)
limit_ZL_ripple  = 0.03;  % peak-to-peak spread of real(ZL) over the band
limit_S11_dB     = -25;   % a uniform line between the measurement planes
limit_S21_dB     = 0.3;   % lossless line, |S21| = 1

% geometry, in mm
unit         = 1e-3;
MSL_length   = 100;
MSL_width    = 2;
MSL_height   = 1;        % strip height over the ground plane
air_height   = 10;
air_width    = 15;

f_max = 5e9;
f     = linspace(1e9, f_max, 201);

Sim_Path = ts_sim_path(mfilename('fullpath'));
Sim_CSX  = 'msl.xml';

%% FDTD setup
FDTD = InitFDTD('NrTS', 20000, 'EndCriteria', 1e-5);
FDTD = SetGaussExcite(FDTD, f_max/2, f_max/2);
% x: propagation direction, terminated by PML; zmin: PEC ground plane
FDTD = SetBoundaryCond(FDTD, {'PML_8','PML_8','MUR','MUR','PEC','MUR'});

%% mesh
CSX = InitCSX();
resolution = c0/f_max/unit/30;   % lambda/30
mesh.x = SmoothMeshLines([-MSL_length/2 0 MSL_length/2], resolution);
mesh.y = SmoothMeshLines([-air_width -MSL_width/2 0 MSL_width/2 air_width], resolution);
mesh.z = SmoothMeshLines([linspace(0,MSL_height,5) air_height], resolution);
CSX = DefineRectGrid(CSX, unit, mesh);

%% the line, as two transmission line ports facing each other
CSX = AddMetal(CSX, 'PEC');
[CSX,port{1}] = AddMSLPort(CSX, 999, 1, 'PEC', ...
    [-MSL_length/2 -MSL_width/2 MSL_height], [0 MSL_width/2 0], 'x', [0 0 -1], ...
    'ExcitePort', true, 'FeedShift', 10*resolution, 'MeasPlaneShift', MSL_length/3);
[CSX,port{2}] = AddMSLPort(CSX, 999, 2, 'PEC', ...
    [ MSL_length/2 -MSL_width/2 MSL_height], [0 MSL_width/2 0], 'x', [0 0 -1], ...
    'MeasPlaneShift', MSL_length/3);

%% run
WriteOpenEMS([Sim_Path '/' Sim_CSX], FDTD, CSX);
Settings.LogFile = [Sim_Path '/openEMS.log'];
Settings.Silent  = opt.Silent;
RunOpenEMS(Sim_Path, Sim_CSX, opt.openEMS_opts, Settings);

%% analysis
port = calcPort(port, Sim_Path, f);

ZL   = port{1}.ZL;
beta = port{1}.beta;
k0   = 2*pi*f/c0;

S11 = port{1}.uf.ref ./ port{1}.uf.inc;
S21 = port{2}.uf.ref ./ port{1}.uf.inc;

% the actual separation of the two measurement planes: MeasPlaneShift is
% requested in drawing units but snaps to the nearest mesh line
meas_distance = (MSL_length - port{1}.measplanepos - port{2}.measplanepos) * unit;

if opt.Plots
    figure
    plot(f/1e9, real(beta)./k0, 'k-', 'LineWidth', 2);
    xlabel('Frequency (GHz)'); ylabel('\beta / k_0');
    title('normalized propagation constant');

    figure
    plot(f/1e9, [real(ZL); imag(ZL)], 'LineWidth', 2);
    xlabel('Frequency (GHz)'); ylabel('Z_L (Ohm)');
    legend({'real','imag'});
    title('characteristic line impedance');

    figure
    plot(f/1e9, 20*log10(abs([S11; S21])), 'LineWidth', 2);
    xlabel('Frequency (GHz)'); ylabel('S-parameter (dB)');
    legend({'S11','S21'});
end

pass = ts_check_rel('propagation constant beta = k0', real(beta), k0, limit_beta);

pass = ts_check('line impedance is real', ...
                max(abs(imag(ZL))./real(ZL)) <= limit_ZL_imag, ...
                'max. |imag(ZL)|/real(ZL) = %.3g (limit %.3g)', ...
                max(abs(imag(ZL))./real(ZL)), limit_ZL_imag) && pass;

ripple = (max(real(ZL)) - min(real(ZL))) / mean(real(ZL));
pass = ts_check('line impedance is frequency independent', ...
                ripple <= limit_ZL_ripple, ...
                'ZL = %.1f .. %.1f Ohm, spread %.2g (limit %.2g)', ...
                min(real(ZL)), max(real(ZL)), ripple, limit_ZL_ripple) && pass;

pass = ts_check('no reflection on the uniform line', ...
                max(20*log10(abs(S11))) <= limit_S11_dB, ...
                'max. S11 = %.1f dB (limit %.1f dB)', ...
                max(20*log10(abs(S11))), limit_S11_dB) && pass;

pass = ts_check('lossless transmission', ...
                max(abs(20*log10(abs(S21)))) <= limit_S21_dB, ...
                'S21 = %.2f .. %.2f dB (limit +/-%.2f dB)', ...
                min(20*log10(abs(S21))), max(20*log10(abs(S21))), limit_S21_dB) && pass;

% the phase of S21 has to match the distance between the measurement planes;
% the measured phase is lifted onto the branch of the reference
phase_ref  = -k0 * meas_distance;
phase_meas = unwrap(angle(S21));
phase_meas = phase_meas + 2*pi*round((phase_ref(1) - phase_meas(1))/(2*pi));
pass = ts_check_rel('S21 phase = -k0 * measurement plane distance', ...
                    phase_meas, phase_ref, limit_beta) && pass;

pass = ts_finish(opt, mfilename, pass, Sim_Path);
