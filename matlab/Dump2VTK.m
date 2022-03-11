function Dump2VTK(filename, fields, mesh, fieldname, varargin)
% Dump2VTK(filename, fields, mesh, fieldname, varargin)
%
%   Dump fields extracted from an hdf5 file to a vtk file format
%
%   possible arguments:
%       'NativeDump': 0 (default) / 1, dump in native coordinate system
%       'CloseAlpha': 0 (default) / 1, repeat first/last line in
%                     alpha-direction for a full cylindrical mesh
%
%   example:
%
%
% openEMS matlab interface
% -----------------------
% author: Thorsten Liebig
%
% See also ReadHDF5FieldData ReadHDF5Mesh GetField_TD2FD GetField_Interpolation

NativeDump = 0;
CloseAlpha = 0;

for n=1:2:numel(varargin)
	if (strcmp(varargin{n},'NativeDump')==1);
		NativeDump =  varargin{n+1};
    elseif (strcmp(varargin{n},'CloseAlpha')==1);
		CloseAlpha =  varargin{n+1};
    end
end

x = mesh.lines{1};
y = mesh.lines{2};
z = mesh.lines{3};

fid = fopen(filename,'w+');

% set nan values to zero
ind = find(isnan(fields));
if (~isempty(ind))
    warning('openEMS:Dump2VTK','field contains nan, setting to zero');
    fields(ind)=0;
end

% set inf values to zero
ind = find(isinf(fields));
if (~isempty(ind))
    warning('openEMS:Dump2VTK','field contains inf, setting to zero');
    fields(ind)=0;
end

if ((CloseAlpha~=0) && (mesh.type==1) && (range(y)<2*pi))
    y(end+1) = y(1)+2*pi;
    fields(:,end+1,:,:) = fields(:,1,:,:);
end

if (mesh.type==0) %write cartesian mesh to vtk
    fprintf(fid,'# vtk DataFile Version 2.0\n');
    fprintf(fid,'Rectilinear Grid by matlab-interface of openEMS\n');
    fprintf(fid,'ASCII\n');
    fprintf(fid,'DATASET RECTILINEAR_GRID\n');

    fprintf(fid,'DIMENSIONS %d %d %d\n',numel(x),numel(y),numel(z));

    fprintf(fid,'X_COORDINATES %d double\n',numel(x));
    fprintf(fid,'%e',x(1));
    for n=2:numel(x)
        fprintf(fid,' %e',x(n));
    end
    fprintf(fid,'\n');

    fprintf(fid,'Y_COORDINATES %d double\n',numel(y));
    fprintf(fid,'%e',y(1));
    for n=2:numel(y)
        fprintf(fid,' %e',y(n));
    end
    fprintf(fid,'\n');

    fprintf(fid,'Z_COORDINATES %d double\n',numel(z));
    fprintf(fid,'%e',z(1));
    for n=2:numel(z)
        fprintf(fid,' %e',z(n));
    end

elseif (mesh.type==1) %write cylindrical mesh to vtk
    fprintf(fid,'# vtk DataFile Version 3.0\n');
    fprintf(fid,'Structured Grid by matlab-interface of openEMS\n');
    fprintf(fid,'ASCII\n');
    fprintf(fid,'DATASET STRUCTURED_GRID\n');

    fprintf(fid,'DIMENSIONS %d %d %d\n',numel(x),numel(y),numel(z));

    fprintf(fid,'POINTS %d double\n',numel(x)*numel(y)*numel(z));
    
    for nz=1:numel(z)
        for ny=1:numel(y)
            for nx=1:numel(x)
                fprintf(fid,'%e %e %e\n',x(nx)*cos(y(ny)),x(nx)*sin(y(ny)),z(nz));
            end
        end
    end
    if ((ndims(fields)==4) && (NativeDump==0))
        [R A Z] = ndgrid(x,y,z);
        sinA = sin(A);
        cosA = cos(A);
        field_CC(:,:,:,1) = fields(:,:,:,1) .* cosA - fields(:,:,:,2) .* sinA;
        field_CC(:,:,:,2) = fields(:,:,:,1) .* sinA + fields(:,:,:,2) .* cosA;
        field_CC(:,:,:,3) = fields(:,:,:,3);
        fields = field_CC;
        clear R A Z sinA cosA field_CC
    end
elseif (mesh.type==2) %write spherical mesh to vtk
    fprintf(fid,'# vtk DataFile Version 3.0\n');
    fprintf(fid,'Structured Grid by matlab-interface of openEMS\n');
    fprintf(fid,'ASCII\n');
    fprintf(fid,'DATASET STRUCTURED_GRID\n');

    fprintf(fid,'DIMENSIONS %d %d %d\n',numel(x),numel(y),numel(z));

    fprintf(fid,'POINTS %d double\n',numel(x)*numel(y)*numel(z));

    for nz=1:numel(z)
        for ny=1:numel(y)
            for nx=1:numel(x)
                fprintf(fid,'%e %e %e\n',...
                    x(nx)*sin(y(ny))*cos(z(nz)),...
                    x(nx)*sin(y(ny))*sin(z(nz)),...
                    x(nx)*cos(y(ny)));
            end
        end
    end

    if ((ndims(fields)==4) && (NativeDump==0))
        [R T A] = ndgrid(x,y,z);
        sinA = sin(A);
        cosA = cos(A);
        sinT = sin(T);
        cosT = cos(T);
        field_CC(:,:,:,1) = fields(:,:,:,1) .* sinT .* cosA + fields(:,:,:,2) .*cosT .* cosA - fields(:,:,:,3) .* sinA;
        field_CC(:,:,:,2) = fields(:,:,:,1) .* sinT .* cosA + fields(:,:,:,2) .*cosT .* sinA + fields(:,:,:,3) .* cosA;
        field_CC(:,:,:,3) = fields(:,:,:,1) .* cosT - fields(:,:,:,2) .*sinT;
        fields = field_CC;
        clear R A T sinA cosA sinT cosT field_CC
    end
end

    
fprintf(fid,'\n\n');

fprintf(fid,'POINT_DATA %d\n',numel(x)*numel(y)*numel(z));
% dump vector field data
if (size(fields,4)>1)
    if (nargin>3)
        fprintf(fid,['VECTORS ' fieldname ' double\n']);
    else
        fprintf(fid,'VECTORS field double\n');
    end
    fclose(fid);
    field_x = fields(:,:,:,1);
    field_y = fields(:,:,:,2);
    field_z = fields(:,:,:,3);
    clear fields
    dumpField(:,1) = field_x(:);
    dumpField(:,2) = field_y(:);
    dumpField(:,3) = field_z(:);
    save('-ascii','-append',filename,'dumpField')
    return
elseif (size(fields,4)==1) % scalar field
    if (nargin>3)
        fprintf(fid,['SCALARS ' fieldname ' double 1\nLOOKUP_TABLE default\n']);
    else
        fprintf(fid,'SCALARS field double 1\nLOOKUP_TABLE default\n');
    end
    fclose(fid);
    dumpField = fields(:);
    save('-ascii','-append',filename,'dumpField')
    return
end

fclose(fid);
