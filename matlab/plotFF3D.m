function h = plotFF3D(nf2ff,varargin)
%  h = plotFF3D(nf2ff,varargin)
%
%  plot normalized 3D far field pattern
%
% input:
%   nf2ff:      output of CalcNF2FF
%
% variable input:
%   'freq_index':  - use the given frequency index, see nf2ff.freq
%                  - default is 1
%   'logscale':    - if set, show far field with logarithmic scale
%                  - set the dB value for point of origin
%                  - values below will be clamped
%   'normalize':   - true/false, normalize linear plot
%                  - default is false, log-plot is always normalized!
%
%   example:
%       plotFF3D(nf2ff, 'freq_index', 2, 'logscale', -20)
%
%       see examples/antennas/infDipol.m
%
% See also CalcNF2FF, plotFFdB, polarFF
% 
% openEMS matlab interface
% -----------------------
% author: Thorsten Liebig, Stefan Mahr

% defaults
logscale = [];
freq_index = 1;
normalize = 0;

for n=1:2:numel(varargin)
    if (strcmp(varargin{n},'logscale')==1);
        logscale = varargin{n+1};
    elseif (strcmp(varargin{n},'freq_index')==1);
        freq_index = varargin{n+1};
    elseif (strcmp(varargin{n},'normalize')==1);
        normalize = varargin{n+1};
    else
        warning('openEMS:plotFF3D',['unknown argument key: ''' varargin{n} '''']);
    end
end

if ((normalize~=0) || ~isempty(logscale))
    E_far = nf2ff.E_norm{freq_index} / max(nf2ff.E_norm{freq_index}(:));
else
    E_far = nf2ff.E_norm{freq_index};
end;

if ~isempty(logscale)
    E_far = 20*log10(E_far)/-logscale + 1;
    E_far = E_far .* ( E_far > 0 );
    titletext = sprintf('electrical far field [dB] @ f = %e Hz',nf2ff.freq(freq_index));
elseif (normalize==0)
    titletext = sprintf('electrical far field [V/m] @ f = %e Hz',nf2ff.freq(freq_index));
else
    titletext = sprintf('normalized electrical far field @ f = %e Hz',nf2ff.freq(freq_index));
end

[theta,phi] = ndgrid(nf2ff.theta,nf2ff.phi);
x = E_far .* sin(theta) .* cos(phi);
y = E_far .* sin(theta) .* sin(phi);
z = E_far .* cos(theta);
%figure
h = surf( x,y,z, E_far );
set(h,'EdgeColor','none');
axis equal
axis off

try
    if (isOctave && (strcmp(graphics_toolkit,'gnuplot')==1))
        warning('openEMS:plotFF3D','Colorbar doesn''t work properly with octave and gnuplot. On problems, try ''colorbar off''');
    end
end

if ~isempty(logscale)
    colorbar('YTick', linspace(0,max(E_far(:)),9), ...
    'YTickLabel',num2str(linspace(logscale, 10*log10(nf2ff.Dmax(freq_index)),9)'));
else
    colorbar;
end

title( titletext );

if (nargout == 0)
  clear h;
end

end
