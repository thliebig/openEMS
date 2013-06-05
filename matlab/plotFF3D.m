function h = plotFF3D(nf2ff,varargin)
%  h = plotFF3D(nf2ff,varargin)
%
%  plot normalized 3D far field pattern
%
% input:
%   nf2ff:      output of CalcNF2FF
%
% variable input:
%   'cellelement': - use element from cell array
%                  - default is 1
%   'logscale':    - if set, show farfield with logarithmic scale
%                  - set the dB value for point of origin
%                  - values below will be clamped
%
%   example:
%       plotFF3D(nf2ff, 'cellelement', 2, 'logscale', -20)
%
%       see examples/NF2FF/infDipol.m
%
% See also CalcNF2FF
% 
% openEMS matlab interface
% -----------------------
% author: Thorsten Liebig, Stefan Mahr

% defaults
logscale = [];
cellelement = 1;

for n=1:2:numel(varargin)
    if (strcmp(varargin{n},'logscale')==1);
        logscale = varargin{n+1};
    elseif (strcmp(varargin{n},'cellelement')==1);
        cellelement = varargin{n+1};
    end
end

E_far_normalized = nf2ff.E_norm{cellelement} / max(nf2ff.E_norm{cellelement}(:));

if ~isempty(logscale)
    E_far_normalized = 20*log10(E_far_normalized)/-logscale + 1;
    ind = find ( E_far_normalized < 0 );
    E_far_normalized(ind) = 0;
    titletext = sprintf('electrical far field [dB] @ f = %e Hz',nf2ff.freq(cellelement));
else
    titletext = sprintf('electrical far field [V/m] @ f = %e Hz',nf2ff.freq(cellelement));
end

[theta,phi] = ndgrid(nf2ff.theta,nf2ff.phi);
x = E_far_normalized .* sin(theta) .* cos(phi);
y = E_far_normalized .* sin(theta) .* sin(phi);
z = E_far_normalized .* cos(theta);
%figure
h = surf( x,y,z, E_far_normalized );
set(h,'EdgeColor','none');
axis equal

title( titletext );
xlabel( 'x' );
ylabel( 'y' );
zlabel( 'z' );

if (nargout == 0)
  clear h;
end

end
