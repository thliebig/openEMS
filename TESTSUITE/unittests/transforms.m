function pass = transforms(varargin)
%pass = transforms(<key>, <value>, ...)
%
% FFT_time2freq and DFT_time2freq against the analytic spectrum of a Gaussian
% pulse and of a sine wave. Every test in this suite that looks at a spectrum
% goes through one of them, so a change of their scaling or of their
% single-sided convention has to be caught here and not through a shifted
% resonance three tests later.
%
% No simulation, so this runs in well under a second.
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

% tolerances: both transforms are exact up to floating point round-off, the
% comparison against the analytic spectrum is limited by the truncation of the
% sampled pulse
tol_analytic = 1e-6;
tol_exact    = 1e-10;

% Gaussian pulse, sampled so that it has decayed at both ends of the record
N   = 2048;            % power of two: FFT_time2freq does not zero-pad
dt  = 1e-12;
tau = 100e-12;
t0  = 1000e-12;        % the pulse has to have decayed at both ends
t   = (0:N-1)*dt;
val = exp(-((t-t0)/tau).^2);

% single-sided analytic spectrum of that pulse (the factor 2 is the
% single-sided convention of both transforms)
X = @(f) 2 * tau*sqrt(pi) * exp(-(pi*f*tau).^2) .* exp(-1j*2*pi*f*t0);

f = linspace(0, 4e9, 51);
pass = ts_check_rel('DFT_time2freq, Gaussian pulse', ...
                    DFT_time2freq(t, val, f), X(f), tol_analytic);

[f_fft, val_fft] = FFT_time2freq(t, val);
% only compare where the spectrum still carries energy
idx = f_fft <= 4e9;
pass = ts_check_rel('FFT_time2freq, Gaussian pulse', ...
                    val_fft(idx), X(f_fft(idx)), tol_analytic) && pass;

% both transforms compute the same sum, so they have to agree exactly
pass = ts_check_rel('FFT_time2freq == DFT_time2freq', ...
                    val_fft(idx), DFT_time2freq(t, val, f_fft(idx)), tol_exact) && pass;

% dropping the leading, still quiet part of a record must not change the
% spectrum -- that is what the phase correction for t(1) is there for, and
% what every test that cuts the excitation out of a time series relies on
i0 = 501;
[f_cut, val_cut] = FFT_time2freq(t(i0:end), val(i0:end));
idx_cut = f_cut <= 4e9;
pass = ts_check_rel('FFT_time2freq, record not starting at t=0', ...
                    val_cut(idx_cut), X(f_cut(idx_cut)), tol_analytic) && pass;

% periodic mode: an integer number of periods of a sine recovers its amplitude
f0  = 1e9;
amp = 0.9;
t   = (0:999)*dt;              % exactly one period of f0
val = amp * sin(2*pi*f0*t);
pass = ts_check_rel('DFT_time2freq, periodic sine amplitude', ...
                    abs(DFT_time2freq(t, val, f0, 'periodic')), amp, tol_exact) && pass;

if opt.Plots
    figure
    plot(f_fft(idx)/1e9, abs(val_fft(idx)), 'k-', f/1e9, abs(X(f)), 'r--', 'LineWidth', 2);
    xlabel('Frequency (GHz)'); ylabel('|X(f)|');
    legend({'FFT\_time2freq','analytic'});
    title('Gaussian pulse spectrum');
end

pass = ts_finish(opt, mfilename, pass);
