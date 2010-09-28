function Dump2VTK(filename, fields, mesh, fieldname)
% Dump2VTK(filename, fields, mesh, fieldname)
%
%   Dump fields extraced from an hdf5 file to a vtk file format
%
%   example:
%
%
% openEMS matlab interface
% -----------------------
% author: Thorsten Liebig
%
% See also ReadHDF5FieldData ReadHDF5Mesh GetField_TD2FD GetField_Interpolation

x = mesh.lines{1};
y = mesh.lines{2};
z = mesh.lines{3};

fid = fopen(filename,'w+');
    
if (mesh.type==0) %write cartesian mesh to vtk
    fprintf(fid,'# vtk DataFile Version 2.0\n');
    fprintf(fid,'Rectilinear Grid by matlab-interface of openEMS\n');
    fprintf(fid,'ASCII\n');
    fprintf(fid,'DATASET RECTILINEAR_GRID\n');

    fprintf(fid,'DIMENSIONS %d %d %d\n',numel(x),numel(y),numel(z));

    fprintf(fid,'X_COORDINATES %d float\n',numel(x));
    fprintf(fid,'%f',x(1));
    for n=2:numel(x)
        fprintf(fid,' %f',x(n));
    end
    fprintf(fid,'\n');

    fprintf(fid,'Y_COORDINATES %d float\n',numel(y));
    fprintf(fid,'%f',y(1));
    for n=2:numel(y)
        fprintf(fid,' %f',y(n));
    end
    fprintf(fid,'\n');

    fprintf(fid,'Z_COORDINATES %d float\n',numel(z));
    fprintf(fid,'%f',z(1));
    for n=2:numel(z)
        fprintf(fid,' %f',z(n));
    end

elseif (mesh.type==1) %write cylindrical mesh to vtk
    fprintf(fid,'# vtk DataFile Version 3.0\n');
    fprintf(fid,'Structured Grid by matlab-interface of openEMS\n');
    fprintf(fid,'ASCII\n');
    fprintf(fid,'DATASET STRUCTURED_GRID\n');

    fprintf(fid,'DIMENSIONS %d %d %d\n',numel(x),numel(y),numel(z));

    fprintf(fid,'POINTS %d float\n',numel(x)*numel(y)*numel(z));
    
    for nz=1:numel(z)
        for ny=1:numel(y)
            for nx=1:numel(x)
                fprintf(fid,'%f %f %f\n',x(nx)*cos(y(ny)),x(nx)*sin(y(ny)),z(nz));
            end
        end
    end
    [R A Z] = ndgrid(x,y,z);
    sinA = sin(A);
    cosA = cos(A);
    field_CC(:,:,:,1) = fields(:,:,:,1) .* cosA - fields(:,:,:,2) .* sinA;
    field_CC(:,:,:,2) = fields(:,:,:,1) .* sinA + fields(:,:,:,2) .* cosA;
    field_CC(:,:,:,3) = fields(:,:,:,3);
    fields = field_CC;
    clear R A Z sinA cosA field_CC
end

    
fprintf(fid,'\n\n');

fprintf(fid,'POINT_DATA %d\n',numel(x)*numel(y)*numel(z));
if (nargin>3)
    fprintf(fid,['VECTORS ' fieldname ' float\n']);
else
    fprintf(fid,'VECTORS field float\n');
end

for nz=1:numel(z)
    for ny=1:numel(y)
        for nx=1:numel(x)
            fprintf(fid,'%f %f %f\n',fields(nx,ny,nz,1),fields(nx,ny,nz,2),fields(nx,ny,nz,3));
        end
    end
end

fclose(fid);
