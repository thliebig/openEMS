function h = plotRefl(port, varargin) 

%  h = plotRefl(port,varargin)
%
%  plot the reflection coefficient of a port into a Smith chart.
%  left and right facing triangles mark the lower and upper cutoff
%  frequency of the pass bands. An asterisk marks the frequnecy with 
%  the lowest reflection.
%
% input:
%   port:      port data structure. Call calcPort with an appropriate 
%              frequency vector before calling this routine
% 
% output:      graphics handle for further modification of the plot.
%
% variable input:
%   'precision':   - number of decimal places (floating point precision) 
%                    for the frequency (always in MHz), default is 0
%   'threshold':   - Threshold value (in dB) for the upper and lower
%                    cutoff frequency, default is -3
%   'fmarkers':    - set lower and upper frequency marker in Hz manually, 
%                    like so: [4e9, 6.5e9]
%   example:
%       myport = calcPort(myport, Sim_Path, linspace(f_0-f_c, f_0+f_c, 200));
%       plotRefl(myport, 'fmarkers', [4e9, 6.5e9]);
%
% See also calcPort
% 
% openEMS matlab interface
% -----------------------
% author: Georg Michel

%defaults
precision = 0;
threshold = -3;
fmarkers = [];

for n=1:2:numel(varargin)
    if (strcmp(varargin{n},'precision')==1);
        precision = varargin{n+1};
    elseif (strcmp(varargin{n},'threshold')==1);
        threshold = varargin{n+1};
    elseif (strcmp(varargin{n},'fmarkers')==1);
        fmarkers = varargin{n+1};
	else
        warning('openEMS:polarFF',['unknown argument key: ''' varargin{n} '''']);
    end
end


if ~isfield(port, 'uf')
  error('Cannot plot the reflection coefficient. Please call calcPort first.');
end

s11 = port.uf.ref ./ port.uf.inc;
ffmt = ['%.', num2str(precision), 'f'];



figure; %new figure

plot([-1, 1], [0, 0], 'k');


axis ([-1.15, 1.15, -1.15, 1.15], "square");
axis off;
hold on

ReZ = [.2; .5; 1; 2]; 
ImZ = 1i * [1 2 5 2];
Z = bsxfun(@plus, ReZ, linspace(-ImZ, ImZ, 256));
Gamma = (Z-1)./(Z+1);
plot(Gamma.', 'k');

ReZ = [.5 .5 1 1 2 2 5 5 10 10];
ImZ = 1i * [-.2; .2; -.5; .5; -1; 1; -2; 2; -5; 5];
Z = bsxfun(@plus, linspace(0, ReZ, 256), ImZ);
Gamma = (Z-1)./(Z+1);
plot(Gamma.', 'k');


angle = linspace (0, 2 * pi, 256); ReZ = [0 5 10];
center = ReZ ./ (ReZ + 1);
radius = 1 ./ (ReZ + 1);
plot(bsxfun(@plus, bsxfun(@times, radius, cos(angle.')), center), bsxfun(@times, radius, sin(angle.')), 'k');


% resistance
ReZ = [0.2 0.5 1 2 5 10]; ImZ = zeros (1, length (ReZ));
rho = (ReZ.^2 + ImZ.^2 - 1 + 2i * ImZ) ./ ((ReZ + 1).^2 + ImZ.^2);

xoffset = [0.1 0.1 0.05 0.05 0.05 0.075]; 
yoffset = -0.03;

for idx = 1:length (ReZ)
    text (real (rho(idx)) - xoffset(idx), ...
        imag (rho(idx)) - yoffset, num2str (ReZ(idx)));
end

% reactance
ReZ = [-0.06 -0.06 -0.06 -0.12 -0.5];
ImZ = [0.2 0.5 1 2 5];


rho = (ReZ.^2 + ImZ.^2 - 1 + 2i * ImZ) ./ ((ReZ + 1).^2 + ImZ.^2);

for idx = 1:length (ImZ)
    text (real (rho(idx)), imag (rho(idx)), [num2str(ImZ(idx)), "j"]);
    text (real (rho(idx)), -imag (rho(idx)), [num2str(-ImZ(idx)), "j"]); end

% zero
rho = (-0.05.^2 + 0.^2 - 1) ./ ((-0.05 + 1).^2 + 0.^2);

text (real (rho), imag (rho), '0');

s11dB = 20*log10(abs(s11));

if(isempty(fmarkers))
upperind = s11dB(1:end-1) < threshold & s11dB(2:end) > threshold;
lowerind = s11dB(1:end-1) > threshold & s11dB(2:end) < threshold;
else
upperind = [nthargout(2, @min, abs(fmarkers(2)-port.f))];
lowerind = [nthargout(2, @min, abs(fmarkers(1)-port.f))];
end

minind = nthargout(2, @min, s11dB);
handle1 = plot(s11(lowerind),['<','b']);
handle2 = plot(s11(upperind),['>','b']);
handle3 = plot(s11(minind),['*', 'b']);
llegend = num2str(port.f(lowerind)(1)/1e6, ffmt);
ulegend = num2str(port.f(upperind)(1)/1e6, ffmt);


if nnz(lowerind) > 1
  for i= 2:nnz(lowerind)
    llegend = strjoin({llegend, num2str(port.f(lowerind)(i)/1e6, ffmt)}, ', ');
  end
end

if nnz(upperind) > 1
  for i= 2:nnz(upperind)
    ulegend = strjoin({ulegend, num2str(port.f(upperind)(i)/1e6, ffmt)}, ', ');
  end
end

legend([handle1, handle2, handle3], {[llegend, " MHz"], ...
                                     [ulegend, " MHz"], ...
                                     [num2str(20*log10(abs(s11(minind))), "%4.0f"), ...
                                    "dB @ ", num2str(port.f(minind)/1e6, ffmt), " MHz"]});
h = plot(s11);

if (nargout == 0)
  clear h;
end

end
