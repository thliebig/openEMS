function nf2ff = ReadNF2FF(nf2ff)
% function nf2ff = ReadNF2FF(nf2ff)
%
% internal function to read calculated nf2ff data, use CalcNF2FF to read
% existing nf2ff data
%
% See also: CalcNF2FF, CreateNF2FFBox
%
% openEMS matlab interface
% -----------------------
% author: Thorsten Liebig, 2012

file = nf2ff.hdf5;

hdf_mesh = ReadHDF5Mesh(file);

nf2ff.r = double(hdf_mesh.lines{1});
nf2ff.theta = double(hdf_mesh.lines{2});
nf2ff.phi = double(hdf_mesh.lines{3});

if isOctave
    nf2ff.freq = double(h5readatt_octave(file,'/nf2ff','Frequency'));
    nf2ff.Prad = double(h5readatt_octave(file,'/nf2ff','Prad'));
    nf2ff.Dmax = double(h5readatt_octave(file,'/nf2ff','Dmax'));
    hdf = load( '-hdf5', file );
    for n=1:numel(nf2ff.freq)
          nf2ff.E_theta{n} = double(hdf.nf2ff.E_theta.FD.(['f' int2str(n-1) '_real']) +1i*hdf.nf2ff.E_theta.FD.(['f' int2str(n-1) '_imag']) );
          nf2ff.E_phi{n} = double(hdf.nf2ff.E_phi.FD.(['f' int2str(n-1) '_real']) +1i*hdf.nf2ff.E_phi.FD.(['f' int2str(n-1) '_imag']) );
          nf2ff.E_norm{n} = double(sqrt(abs(nf2ff.E_theta{n}).^2+abs(nf2ff.E_phi{n}).^2));
    end
else
    nf2ff.freq = double(h5readatt(file,'/nf2ff','Frequency'));
    nf2ff.Prad = double(h5readatt(file,'/nf2ff','Prad'));
    nf2ff.Dmax = double(h5readatt(file,'/nf2ff','Dmax'));

    for n=1:numel(nf2ff.freq)
        nf2ff.E_theta{n} = double(h5read(file,['/nf2ff/E_theta/FD/f' int2str(n-1) '_real']) + 1i*h5read(file,['/nf2ff/E_theta/FD/f' int2str(n-1) '_imag']));
        nf2ff.E_phi{n} = double(h5read(file,['/nf2ff/E_phi/FD/f' int2str(n-1) '_real']) + 1i*h5read(file,['/nf2ff/E_phi/FD/f' int2str(n-1) '_imag']));
        nf2ff.E_norm{n} = double(sqrt(abs(nf2ff.E_theta{n}).^2+abs(nf2ff.E_phi{n}).^2));
    end
end
