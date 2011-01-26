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
    fprintf(fid,'%e',x(1));
    for n=2:numel(x)
        fprintf(fid,' %e',x(n));
    end
    fprintf(fid,'\n');

    fprintf(fid,'Y_COORDINATES %d float\n',numel(y));
    fprintf(fid,'%e',y(1));
    for n=2:numel(y)
        fprintf(fid,' %e',y(n));
    end
    fprintf(fid,'\n');

    fprintf(fid,'Z_COORDINATES %d float\n',numel(z));
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

    fprintf(fid,'POINTS %d float\n',numel(x)*numel(y)*numel(z));
    
    for nz=1:numel(z)
        for ny=1:numel(y)
            for nx=1:numel(x)
                fprintf(fid,'%e %e %e\n',x(nx)*cos(y(ny)),x(nx)*sin(y(ny)),z(nz));
            end
        end
    end
    if (ndims(fields)==4)
        [R A Z] = ndgrid(x,y,z);
        sinA = sin(A);
        cosA = cos(A);
        field_CC(:,:,:,1) = fields(:,:,:,1) .* cosA - fields(:,:,:,2) .* sinA;
        field_CC(:,:,:,2) = fields(:,:,:,1) .* sinA + fields(:,:,:,2) .* cosA;
        field_CC(:,:,:,3) = fields(:,:,:,3);
        fields = field_CC;
        clear R A Z sinA cosA field_CC
    end
end

    
fprintf(fid,'\n\n');

fprintf(fid,'POINT_DATA %d\n',numel(x)*numel(y)*numel(z));
% dump vector field data
if (ndims(fields)==4)
    if (nargin>3)
        fprintf(fid,['VECTORS ' fieldname ' float\n']);
    else
        fprintf(fid,'VECTORS field float\n');
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
end

% dump scalar field data
if (ndims(fields)==3)
    if (nargin>3)
        fprintf(fid,['SCALARS ' fieldname ' float 1\nLOOKUP_TABLE default\n']);
    else
        fprintf(fid,'SCALARS field float 1\nLOOKUP_TABLE default\n');
    end
    fclose(fid);
    dumpField = fields(:);
    save('-ascii','-append',filename,'dumpField')
    return
end

fclose(fid);
