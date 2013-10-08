function h = polarFF(nf2ff,varargin)
%  h = polarFF(nf2ff,varargin)
%
%  plot polar far field pattern
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
%   'normalize':   - true/false, normalize linear plot
%                  - default is false, log-plot is always normalized!
%   'logscale':    - if set, plot logarithmic polar
%                  - set the dB value for point of origin if scalar
%                  - set point of origin and maximum if 2-element array
%                  - values below minimum will be clamped
%                  - default is -20
%   'xtics':       - set the number of tics for polar grid
%                  - default is 5
%
%   example:
%       polarFF(nf2ff, 'freq_index', 2, ...
%                      'xaxis', 'phi', 'param', [1 46 91] );
%
%       polarFF(..., 'normalize', true );
%       polarFF(..., 'logscale', -30 );
%       polarFF(..., 'logscale', [-30 10]);
%
%       polarFF(..., 'xtics', 10);
%
%       see examples/antenna/infDipol.m
%
% See also CalcNF2FF, plotFFdB, plotFF3D
% 
% openEMS matlab interface
% -----------------------
% author: Thorsten Liebig, Stefan Mahr

% defaults
freq_index = 1;
xaxis = 'phi';
param = 1;
logscale = [];
xtics = 5;
normalize = 0;

for n=1:2:numel(varargin)
    if (strcmp(varargin{n},'freq_index')==1);
        freq_index = varargin{n+1};
    elseif (strcmp(varargin{n},'xaxis')==1);
        xaxis = varargin{n+1};
    elseif (strcmp(varargin{n},'param')==1);
        param = varargin{n+1};
    elseif (strcmp(varargin{n},'normalize')==1);
        normalize = varargin{n+1};
    elseif (strcmp(varargin{n},'logscale')==1);
        logscale = varargin{n+1};
    elseif (strcmp(varargin{n},'xtics')==1);
        xtics = varargin{n+1};
    else
        warning('openEMS:polarFF',['unknown argument key: ''' varargin{n} '''']);
    end
end

E_far_max = max(nf2ff.E_norm{freq_index}(:));
if ~isempty(logscale)
    gridmin = logscale(1);

    Dmax = 10*log10(nf2ff.Dmax(freq_index));
    E_far_scale = Dmax - gridmin;
    E_far = 20*log10(nf2ff.E_norm{freq_index}) - 20*log10(E_far_max) + E_far_scale;
    E_far = E_far .* ( E_far > 0 );
    E_far = E_far ./ E_far_scale;

    titletext = sprintf('electrical far field [dBi] @ f = %e Hz',nf2ff.freq(freq_index));

    if numel(logscale) == 2 % normalize to maximum grid
        gridmax = logscale(2);
        E_far = E_far .* E_far_scale/(gridmax-gridmin);
    else
        gridmax = Dmax;
    end
elseif (normalize==0)
    E_far = nf2ff.E_norm{freq_index};

    titletext = sprintf('electrical far field [V/m] @ f = %e Hz',nf2ff.freq(freq_index));

    gridmin = 0;
    gridmax = E_far_max;
else % normalize == 1
    E_far = nf2ff.E_norm{freq_index} / E_far_max;

    titletext = sprintf('normalized electrical far field @ f = %e Hz',nf2ff.freq(freq_index));

    gridmin = 0;
    gridmax = 1;
end


if (strcmp(xaxis,'theta')==1);
    xax = nf2ff.theta(:);
    yax = E_far(:,param);
    parval = nf2ff.phi(param);
    param = 'phi';
elseif (strcmp(xaxis,'phi')==1);
    xax = nf2ff.phi(:);
    yax = E_far(param,:)';
    parval = nf2ff.theta(param);
    param = 'theta';
else
    error('openEMS:polarFF','unknown parameter to ''xaxis''');
end

if ~isempty(logscale)
    scalegrid = 1;
else
    scalegrid = gridmax;
end

% workaround for polar plot
gridcolor = [0.85 0.85 0.85];
% plot xtics circles
a=linspace(0,2*pi,60);
b=linspace(0,scalegrid,xtics+1);
b=repmat(b(2:end),numel(a),1)';
a=repmat(a,size(b,1),1);
[x,y] = pol2cart(a,b);
h = plot(x',y');
%h=polar(a,b,'-k');
set(h,'Color',gridcolor);
set(h(end),'Color',gridcolor*0.8);
hold on;
% plot degree lines
a=bsxfun(@plus,[0:pi/6:pi-pi/6],[0 pi]');
b=scalegrid.*ones(size(a));
h=polar(a,b,'-k');
set(h,'Color',gridcolor);
set(h([1 4]),'Color',gridcolor*0.8);
text(scalegrid*0.05,scalegrid*0.05,num2str(gridmin))
text(scalegrid*1.05,scalegrid*0.05,num2str(gridmax))


% draw far field
xax = repmat(xax,1,size(yax,2));
[x,y] = pol2cart(xax,yax);
h = plot(x,y);
%h = polar( xax, yax );

% legend
ylabel( sprintf('%s / deg', xaxis) );
title( titletext );
createlegend = @(d)sprintf('%s = %3.1f',param,d / pi * 180);
legendtext = arrayfun(createlegend,parval,'UniformOutput',0);
legend( h, legendtext ,'location','southeast');

% workaround for polar plot
axis equal tight
axis ([-scalegrid scalegrid -scalegrid scalegrid]);
axis off
hold off
set(gcf,'Color','white');

if (nargout == 0)
  clear h;
end

end
