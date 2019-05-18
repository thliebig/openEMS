function [delay, fidelity, nf2ff_out] = DelayFidelity(nf2ff, port, path, weight_theta, weight_phi, theta, phi, f_0, f_c, varargin)
% [delay, fidelity] = DelayFidelity(nf2ff, port, path, weight_theta, weight_phi, theta, phi, f_lo, f_hi, varargin)
% 
% 
% This function calculates the time delay from the source port to the phase center of the antenna and the fidelity. 
% The fidelity is the similarity between the excitation pulse and the radiated pulse (normalized scalar product).
% The resolution of the delay will be equal to or better than ((f_0 + f_c)*Oversampling)^-1 when using Gaussian excitation.
% Oversampling is an input parameter to InitFDTD. The rows of delay and fidelity correspond to theta and the columns to phi.
% 
% input:
%   nf2ff: return value of CreateNF2FFBox.
%   port: return value of AddLumpedPort
%   path: path of the simulation results.
%   weight_theta: weight of the E_theta component 
%   weight_phi: weight of the E_phi component
%     -> with both (possibly complex) parameters any polarization can be examined
%   theta: theta values to be simulated
%   phi: phi values to be simulated
%   f_0: center frequency of SetGaussExcite
%   f_c: cutoff frequency of SetGaussExcite
%
% variable input:
%   'Center': phase center of the antenna for CalcNF2FF
%   'Radius': radius for CalcNF2FF
%   'Mode': mode CalcNF2FF
%   
% example:
%   theta = [-180:10:180] * pi / 180;
%   phi = [0, 90] * pi / 180;
%   % use circular right handed polarization
%   [delay, fidelity] = DelayFidelity2(nf2ff, port, Sim_Path, -1i, 1, theta, phi, f_0, f_c, 'Mode', 1);
%   figure
%   polar(theta.', delay(:,1) * 3e11); % delay in mm
%   figure
%   polar(theta', (fidelity(:,1)-0.95)/0.05); % last 5 percent of fidelity
% 
% Author: Georg Michel 

C0 = 299792458;
center = [0, 0, 0];
radius = 1;
nf2ff_mode = 0;

for n=1:2:numel(varargin)
    if (strcmp(varargin{n},'Center')==1);
        center = varargin{n+1};
    elseif (strcmp(varargin{n},'Radius')==1);
        radius = varargin{n+1};
    elseif (strcmp(varargin{n},'Mode')==1);
        nf2ff_mode = varargin{n+1};
    end
end


port_ut = load(fullfile(path, port.U_filename));
port_it = load(fullfile(path, port.I_filename));
dt = port_ut(2,1) - port_ut(1,1);
fftsize = 2^(nextpow2(size(port_ut)(1)) + 1);
df = 1 / (dt * fftsize);
uport = fft(port_ut(:, 2), fftsize)(1:fftsize/2+1);
iport = fft(port_it(:, 2), fftsize)(1:fftsize/2+1);
fport = df * (0:fftsize/2);
f_ind = find(fport > (f_0 - f_c ) & fport < (f_0 + f_c));
disp(["frequencies: ", num2str(numel(f_ind))]);
exc_f = uport.' + iport.' * port.Feed_R; %excitation in freq domain
exc_f(!f_ind) = 0;
exc_f /= sqrt(exc_f * exc_f'); % normalization (transposing also conjugates)

nf2ff = CalcNF2FF(nf2ff, path, fport(f_ind), theta, phi, ...
        'Center', center, 'Radius', radius, 'Mode', nf2ff_mode);
radfield = weight_theta * cell2mat(nf2ff.E_theta) + weight_phi * cell2mat(nf2ff.E_phi); % rows: theta(f1), columns: phi(f1), phi(f2), ...phi(fn)
radfield = reshape(radfield, [length(nf2ff.theta), length(nf2ff.phi), length(nf2ff.freq)]);
correction = reshape(exp(-2i*pi*nf2ff.r/C0*nf2ff.freq), 1,1,numel(nf2ff.freq)); %dimensions: theta, phi, frequencies
radfield = radfield./correction; % correct for radius delay
% normalize radfield
radnorm = sqrt(dot(radfield, radfield, 3));
radfield ./= radnorm;

%initialize radiated field in fully populated frequency domain
rad_f = zeros([numel(nf2ff.theta), numel(nf2ff.phi), numel(fport)]);
rad_f(:, :, f_ind) = radfield; % assign selected frequencies
exc_f = reshape(exc_f, [1,1,numel(exc_f)]); %make exc_f conformant with rad_f

cr_f = rad_f .* conj(exc_f); % calculate cross correlation
% calculate the cross correlation in time domain (analytic signal)
cr = ifft(cr_f(:, :, 1:end-1), [], 3) * (numel(fport) -1); % twice the FFT normalization (sqrt^2) because product of two normalized functions 
%search for the maximum of the envelope
[fidelity, delay_ind] = max(abs(cr), [], 3);
delay = (delay_ind - 1) * dt * 2; % double time step because of single-sided FFT
nf2ff_out = nf2ff; %possibly needed for plotting the far field and other things
disp(["DelayFidelity: delay resolution = ", num2str(dt*2e9), "ns"]);
return;


