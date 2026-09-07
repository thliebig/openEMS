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

% read attributes
nf2ff.freq = ReadHDF5Attribute(file,'/nf2ff','Frequency');
nf2ff.Prad = ReadHDF5Attribute(file,'/nf2ff','Prad');
nf2ff.Dmax = ReadHDF5Attribute(file,'/nf2ff','Dmax');

try
    nf2ff.Eps_r = ReadHDF5Attribute(file,'/nf2ff','Eps_r');
catch
    nf2ff.Eps_r = ones(size(nf2ff.freq));
end
try
    nf2ff.Mue_r = ReadHDF5Attribute(file,'/nf2ff','Mue_r');
catch
    nf2ff.Mue_r = ones(size(nf2ff.freq));
end

% CalcNF2FF asks the nf2ff binary for the legacy format: the complex data
% split into a real and an imaginary dataset, stored (phi,theta). A file
% written without that request (by the python interface or by calling the
% nf2ff binary directly) holds one compound complex dataset per frequency,
% stored (theta,phi). Octave is unable to read the compound type at all.

if isOctave
    hdf = load( '-hdf5', file );
    if ~isfield(hdf.nf2ff.E_theta.FD, 'f0_real')
        error('openEMS:ReadNF2FF', ['"' file '" holds the far field as a ' ...
              'compound complex type, which Octave cannot read. Delete the ' ...
              'file and recalculate it with CalcNF2FF.']);
    end
    for n=1:numel(nf2ff.freq)
        nf2ff.E_theta{n} = double(hdf.nf2ff.E_theta.FD.(['f' int2str(n-1) '_real']) +1i*hdf.nf2ff.E_theta.FD.(['f' int2str(n-1) '_imag']) );
        nf2ff.E_phi{n} = double(hdf.nf2ff.E_phi.FD.(['f' int2str(n-1) '_real']) +1i*hdf.nf2ff.E_phi.FD.(['f' int2str(n-1) '_imag']) );
        nf2ff.E_norm{n} = double(sqrt(abs(nf2ff.E_theta{n}).^2+abs(nf2ff.E_phi{n}).^2));
        nf2ff.P_rad{n} = double(hdf.nf2ff.P_rad.FD.(['f' int2str(n-1)]));
    end
else
    % matlab compatibility to older versions
    if verLessThan('matlab','7.12')
        % read data
        for n=1:numel(nf2ff.freq)
            nf2ff.E_theta{n} = double(hdf5read(file,['/nf2ff/E_theta/FD/f' int2str(n-1) '_real']) + 1i*hdf5read(file,['/nf2ff/E_theta/FD/f' int2str(n-1) '_imag']));
            nf2ff.E_phi{n} = double(hdf5read(file,['/nf2ff/E_phi/FD/f' int2str(n-1) '_real']) + 1i*hdf5read(file,['/nf2ff/E_phi/FD/f' int2str(n-1) '_imag']));
            nf2ff.E_norm{n} = double(sqrt(abs(nf2ff.E_theta{n}).^2+abs(nf2ff.E_phi{n}).^2));
            nf2ff.P_rad{n} = double(hdf5read(file,['/nf2ff/P_rad/FD/f' int2str(n-1)]));
        end
    else
        info = h5info(file,'/nf2ff/E_theta/FD');
        legacy = any(strcmp({info.Datasets.Name}, 'f0_real'));
        % read data
        for n=1:numel(nf2ff.freq)
            if legacy
                nf2ff.E_theta{n} = double(h5read(file,['/nf2ff/E_theta/FD/f' int2str(n-1) '_real']) + 1i*h5read(file,['/nf2ff/E_theta/FD/f' int2str(n-1) '_imag']));
                nf2ff.E_phi{n} = double(h5read(file,['/nf2ff/E_phi/FD/f' int2str(n-1) '_real']) + 1i*h5read(file,['/nf2ff/E_phi/FD/f' int2str(n-1) '_imag']));
                nf2ff.P_rad{n} = double(h5read(file,['/nf2ff/P_rad/FD/f' int2str(n-1)]));
            else
                nf2ff.E_theta{n} = ReadComplexFD(file,['/nf2ff/E_theta/FD/f' int2str(n-1)]);
                nf2ff.E_phi{n} = ReadComplexFD(file,['/nf2ff/E_phi/FD/f' int2str(n-1)]);
                nf2ff.P_rad{n} = double(h5read(file,['/nf2ff/P_rad/FD/f' int2str(n-1)])).';
            end
            nf2ff.E_norm{n} = double(sqrt(abs(nf2ff.E_theta{n}).^2+abs(nf2ff.E_phi{n}).^2));
        end
    end
end

% Calculation of right- and left-handed circular polarization
% adopted from
% 2012, Tim Pegg <teepegg@gmail.com>

% cleanup (if exist)
nf2ff.E_cprh = [];
nf2ff.E_cplh = [];

% Setup vectors for converting to LHCP and RHCP polarization senses
[THETHA PHI] = ndgrid(nf2ff.theta,nf2ff.phi);
cosphi = cos(PHI);
sinphi = sin(PHI);

for f=1:numel(nf2ff.freq)
    nf2ff.E_cprh{f} = (cosphi+1i*sinphi) .* (nf2ff.E_theta{f}+1i*nf2ff.E_phi{f})/sqrt(2);
    nf2ff.E_cplh{f} = (cosphi-1i*sinphi) .* (nf2ff.E_theta{f}-1i*nf2ff.E_phi{f})/sqrt(2);
end

end

function data = ReadComplexFD(file, dataset)
% read a compound complex dataset, stored (theta,phi), and return it in the
% [theta_idx, phi_idx] order used by all nf2ff results
raw = h5read(file, dataset);
data = double(raw.r).' + 1i*double(raw.i).';
end
