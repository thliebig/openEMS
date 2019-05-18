function hdf_mesh = ReadHDF5Mesh(file)
% function hdf_mesh = ReadHDF5Mesh(file)
%
%   Get the raw mesh data stored in the hdf5 dump file created by openEMS
%
% returns:
% hdf_mesh.type     (0-> Cartesian, 1-> cylindrical mesh type)
% hdf_mesh.names    (e.g. 'Mesh/y')
% hdf_mesh.lines    (e.g. [0,1,2,3,4])
%
% openEMS matlab interface
% -----------------------
% author: Thorsten Liebig
%
% See also ReadHDF5FieldData

isOctave = exist('OCTAVE_VERSION','builtin') ~= 0;
if isOctave
    hdf_mesh = ReadHDF5Mesh_octave(file);
    return
end

info = hdf5info(file);

for n=1:numel(info.GroupHierarchy.Groups)
    if strcmp(info.GroupHierarchy.Groups(n).Name,'/Mesh')
        for m=1:numel(info.GroupHierarchy.Groups(n).Datasets)
            names{m} = info.GroupHierarchy.Groups(n).Datasets(m).Name;
        end
    end
end

hdf_mesh.names = names;
for n=1:numel(names)
    hdf_mesh.lines{n} = double(hdf5read(file,names{n}));
end

if (strcmp(names{1},'/Mesh/alpha'))
    % alpha and rho are in the wrong order, flip to have rho, alpha, z
    hdf_mesh.names(1:2) = fliplr(hdf_mesh.names(1:2));
    hdf_mesh.lines(1:2) = fliplr(hdf_mesh.lines(1:2));
    hdf_mesh.type=1;
    return
end
if (strcmp(names{1},'/Mesh/phi'))
    % reorder coordinates
    hdf_mesh.names = hdf_mesh.names([2 3 1]);
    hdf_mesh.lines = hdf_mesh.lines([2 3 1]);
    hdf_mesh.type=2;
    return
end

hdf_mesh.type=0;


function hdf_mesh = ReadHDF5Mesh_octave(file)
hdf = load( '-hdf5', file );
hdf_mesh.names = fieldnames(hdf.Mesh);
hdf_mesh.type = 0; % Cartesian mesh
for n=1:numel(hdf_mesh.names)
    hdf_mesh.lines{n} = hdf.Mesh.(hdf_mesh.names{n});
    hdf_mesh.names{n} = ['/Mesh/' hdf_mesh.names{n}];
    if strcmp(hdf_mesh.names{n},'/Mesh/alpha')
        hdf_mesh.type = 1; % cylindrical mesh
    end
    if strcmp(hdf_mesh.names{n},'/Mesh/phi')
        hdf_mesh.type = 2; % cylindrical mesh
    end
end

if (hdf_mesh.type==1)
    % alpha and rho are in the wrong order, flip to have rho, alpha, z
    hdf_mesh.names(1:2) = fliplr(hdf_mesh.names(1:2));
    hdf_mesh.lines(1:2) = fliplr(hdf_mesh.lines(1:2));
end
if (hdf_mesh.type==2)
    % alpha and rho are in the wrong order, flip to have rho, alpha, z
    hdf_mesh.names = hdf_mesh.names([2 3 1]);
    hdf_mesh.lines = hdf_mesh.lines([2 3 1]);
end
