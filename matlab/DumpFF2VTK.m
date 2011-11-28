function DumpFF2VTK(filename, farfield, thetaRange, phiRange, scale)
%  DumpFF2VTK(filename, farfield, thetaRange, phiRange, scale)
%
%  Dump 3D far field pattern to a vtk file
%
%   example:
%       see examples/NF2FF/infDipol.m
%
% See also CreateNF2FFBox, AnalyzeNF2FF
% 
% openEMS matlab interface
% -----------------------
% author: Thorsten Liebig


if (nargin<5)
    scale = 1;
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
dumpField = farfield(:);
save('-ascii','-append',filename,'dumpField')

