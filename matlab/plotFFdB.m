function h = plotFFdB(nf2ff,varargin)
%  h = plotFFdB(nf2ff,varargin)
%
%  plot far field pattern in dBi
%
% input:
%   nf2ff:      output of CalcNF2FF
%
% variable input:
%   'freq_index':  - use the given frequency index, see nf2ff.freq
%                  - default is 1
%   'xaxis':       - 'phi' (default) or 'theta'
%   'param':       - array positions of parametric plot
%                  - if xaxis='phi', theta is parameter, and vice versa
%                  - default is 1
%
%   example:
%       plotFFdB(nf2ff, 'freq_index', 2, ...
%                       'xaxis', 'phi', 'param', [1 46 91])
%
%       see examples/NF2FF/infDipol.m
%
% See also CalcNF2FF, plotFF3D, polarFF
% 
% openEMS matlab interface
% -----------------------
% author: Thorsten Liebig, Stefan Mahr

% defaults
freq_index = 1;
xaxis = 'phi';
param = 1;

for n=1:2:numel(varargin)
    if (strcmp(varargin{n},'freq_index')==1);
        freq_index = varargin{n+1};
    elseif (strcmp(varargin{n},'xaxis')==1);
        xaxis = varargin{n+1};
    elseif (strcmp(varargin{n},'param')==1);
        param = varargin{n+1};
    else
        warning('openEMS:plotFFdB',['unknown argument key: ''' varargin{n} '''']);
    end
end

D_log = nf2ff.E_norm{freq_index} / max(nf2ff.E_norm{freq_index}(:));
D_log = 20*log10(D_log) + 10*log10(nf2ff.Dmax(freq_index));

if (strcmp(xaxis,'theta')==1);
    xax = nf2ff.theta;
    yax = D_log(:,param);
    parval = nf2ff.phi(param);
    param = 'phi';
elseif (strcmp(xaxis,'phi')==1);
    xax = nf2ff.phi;
    yax = D_log(param,:);
    parval = nf2ff.theta(param);
    param = 'theta';
else
    error('openEMS:plotFFdB','unknown parameter to ''xaxis''');
end

%figure
h = plot( xax / pi * 180 , yax );
xlabel( sprintf('%s (deg)',xaxis ));
ylabel( 'directivity (dBi)');

createlegend = @(d)sprintf('%s = %3.1f',param,d / pi * 180);
legendtext = arrayfun(createlegend,parval,'UniformOutput',0);
legend( legendtext );
title( sprintf('far field pattern @ f = %e Hz',nf2ff.freq(freq_index)) );
grid on;

if (nargout == 0)
  clear h;
end

end
