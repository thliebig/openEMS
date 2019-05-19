function DumpFF2VTK(filename, farfield, thetaRange, phiRange, varargin)
%  DumpFF2VTK(filename, farfield, thetaRange, phiRange, varargin)
%
%  Dump 3D far field pattern to a vtk file
%
% input:
%   filename:      filename of VTK file, existing file will be overwritten
%   farfield:      far field in V/m
%   thetaRange:    theta range in deg
%   phiRange:      phi range in deg
%
% variable input:
%   'scale':       - linear scale of plot, doesn't affect gain values
%   'logscale':    - if set, show far field with logarithmic scale
%                  - set the dB value for point of origin
%                  - values below will be clamped
%   'maxgain':     - add max gain in dB to normalized far field
%                  - only valid if logscale is set
%                  - default is 0dB
%
%   example:
%       DumpFF2VTK(filename, farfield, thetaRange, phiRange, ...
%                    'scale', 2, 'logscale', -20, 'maxgain', 3)
%
%       see also examples/NF2FF/infDipol.m
%
% See also CreateNF2FFBox, CalcNF2FF
% 
% openEMS matlab interface
% -----------------------
% author: Thorsten Liebig


% defaults
scale = 1;
maxgain = 0;
logscale = [];

for n=1:2:numel(varargin)
    if (strcmp(varargin{n},'maxgain')==1);
        maxgain = varargin{n+1};
    elseif (strcmp(varargin{n},'logscale')==1);
        logscale = varargin{n+1};
    elseif (strcmp(varargin{n},'scale')==1);
        scale = varargin{n+1};
    end
end

if ~isempty(logscale)
    farfield = 20*log10(farfield) + maxgain - logscale;
    ind = find(farfield<0);
    farfield(ind)=0;
else
    % force 0 for linear plot
    logscale = 0;
end

t = thetaRange*pi/180;
a = phiRange*pi/180;

fid = fopen(filename,'w+');

% set nan values to zero
ind = find(isnan(farfield));
if (~isempty(ind))
    warning('openEMS:Dump2VTK','field contains nan, setting to zero');
    farfield(ind)=0;
end

% set inf values to zero
ind = find(isinf(farfield));
if (~isempty(ind))
    warning('openEMS:Dump2VTK','field contains inf, setting to zero');
    farfield(ind)=0;
end


fprintf(fid,'# vtk DataFile Version 3.0\n');
fprintf(fid,'Structured Grid by matlab-interface of openEMS\n');
fprintf(fid,'ASCII\n');
fprintf(fid,'DATASET STRUCTURED_GRID\n');

fprintf(fid,'DIMENSIONS %d %d %d\n',1,numel(t),numel(a));

fprintf(fid,'POINTS %d double\n',numel(t)*numel(a));

for na=1:numel(phiRange)
    for nt=1:numel(thetaRange)
            fprintf(fid,'%e %e %e\n',...
                scale*farfield(nt,na)*sin(t(nt))*cos(a(na)),...
                scale*farfield(nt,na)*sin(t(nt))*sin(a(na)),...
                scale*farfield(nt,na)*cos(t(nt)));
    end
end



fprintf(fid,'\n\n');

fprintf(fid,'POINT_DATA %d\n',numel(t)*numel(a));

fprintf(fid,['SCALARS gain double 1\nLOOKUP_TABLE default\n']);
fclose(fid);
dumpField = farfield(:) + logscale;
save('-ascii','-append',filename,'dumpField')
