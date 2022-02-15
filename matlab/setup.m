function setup()
% function setup()
%
% setup openEMS Matlab/octave interface
%
% openEMS matlab/octave interface
% -----------------------
% author: Thorsten Liebig (2011-2017)

disp('setting up openEMS matlab/octave interface')

% cd to directory of this file and restore current path at the end
current_path = pwd;
dir  = fileparts( mfilename('fullpath') );
cd(dir);

if isOctave()
    disp('compiling oct files')
    fflush(stdout);
    if isunix
        [res, fn_so] = unix('find /usr/lib     -name libhdf5.so');
        [res, fn_h]  = unix('find /usr/include -name hdf5.h | sort -r | head -1');
        if length(fn_so)>0 && length(fn_h)>0
            [hdf5lib_dir, hdf5lib_fn, ext] = fileparts(fn_so);
            disp(["HDF5 library path found at: " hdf5lib_dir])

            [hdf5inc_dir, hdf5inc_fn, ext] = fileparts(fn_h);
            disp(["HDF5 include path found at: " hdf5inc_dir])
            mkoctfile("h5readatt_octave.cc", ["-L" hdf5lib_dir], ["-I" hdf5inc_dir], "-L hdf5")
        else
            mkoctfile -lhdf5 h5readatt_octave.cc
        end
    else
        mkoctfile -lhdf5 h5readatt_octave.cc
    end
else
    disp('Matlab does not need this function. It is Octave only.')
end

cd(current_path);
