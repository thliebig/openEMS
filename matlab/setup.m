function setup()
% function setup()
%
% setup openEMS Matlab/octave interface
%
% openEMS matlab/octave interface
% -----------------------
% author: Thorsten Liebig (2011)

disp('setting up openEMS matlab/octave interface')

% cd to directory of this file and restore current path at the end
current_path = pwd;
dir  = fileparts( mfilename('fullpath') );
cd(dir);

if isOctave()
    disp('compiling oct files')
    fflush(stdout)
    if isunix
        mkoctfile -lhdf5 -DH5_USE_16_API h5readatt_octave.cc
    else
        mkoctfile -lhdf5 -DH5_USE_16_API h5readatt_octave.cc
    end
else
    disp('Matlab does not need this function. It is Octave only.')
end

cd(current_path);
