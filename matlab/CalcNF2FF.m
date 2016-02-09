function nf2ff = CalcNF2FF(nf2ff, Sim_Path, freq, theta, phi, varargin)
% function nf2ff = CalcNF2FF(nf2ff, Sim_Path, freq, theta, phi, varargin)
%
% Calculate the near-field to far-field transformation created by
% CreateNF2FFBox
%
% IMPORTANT:
% Make sure to define the correct nf2ff phase center, aka. central antenna
% position! See optional parameter below!! Default is [0 0 0]
%
% parameter:
% nf2ff:    data structure created by CreateNF2FFBox
% Sim_Path: path to simulation data
% freq:     array of frequencies to analyse
% theta,phi: spherical coordinates to evaluate the far-field on (in radians)
%
% optional paramater:
% 'Center': nf2ff phase center, default is [0 0 0]
%           !! Make sure the center is never outside of your nf2ff box!!
%           Definition is the correct coordinate system necessary
%           --> either Cartesian or cylindrical coordinates
% 'Mode':   'Mode', 0 -> read only, if data already exist (default)
%           'Mode', 1 -> calculate anyway, overwrite existing
%           'Mode', 2 -> read only, fail if not existing
% 'Outfile': alternative nf2ff result hdf5 file name
%            default is: <nf2ff.name>.h5
% 'Verbose': set verbose level for the nf2ff calculation 0-2 supported
% 'Radius':  specify the radius for the nf2ff
% 'Eps_r':   specify the relative electric permittivity for the nf2ff
% 'Mue_r':   specify the relative magnetic permeability for the nf2ff
%
% 'Mirror':  Add mirroring in a given direction (dir), with a given 
%            mirror type (PEC or PMC) and a mirror position in the given
%            direction.
%            Example: 'Mirror', {0, 'PMC', +100}
%
% See also: CreateNF2FFBox, ReadNF2FF
%
% openEMS matlab interface
% -----------------------
% author: Thorsten Liebig, 2012

mode = 0;

filename = nf2ff.name;
nf2ff_xml.Planes = {};

nf2ff_xml.ATTRIBUTE.Outfile = [filename '.h5'];

if (isfield(nf2ff,'Eps_r'))
    nf2ff_xml.ATTRIBUTE.Eps_r = nf2ff.Eps_r;
end
if (isfield(nf2ff,'Mue_r'))
    nf2ff_xml.ATTRIBUTE.Mue_r = nf2ff.Mue_r;
end

for n=1:2:numel(varargin)-1
    if (strcmp(varargin{n},'Mode'))
        mode = varargin{n+1};
    elseif (strcmp(varargin{n},'Mirror'))
        if isfield(nf2ff_xml,'Mirror')
            pos = length(nf2ff_xml.Mirror)+1;
        else
            pos = 1;
        end
        nf2ff_xml.Mirror{pos}.ATTRIBUTE.Dir=varargin{n+1}{1};
        nf2ff_xml.Mirror{pos}.ATTRIBUTE.Type=varargin{n+1}{2};
        nf2ff_xml.Mirror{pos}.ATTRIBUTE.Pos=varargin{n+1}{3};
    else
        nf2ff_xml.ATTRIBUTE.(varargin{n})=varargin{n+1};
    end
end

for (n=1:numel(nf2ff.filenames_E))
    if (nf2ff.directions(n)~=0)
        files_E = dir([Sim_Path '/*' nf2ff.filenames_E{n} '.h5']);
        files_H = dir([Sim_Path '/*' nf2ff.filenames_H{n} '.h5']);
        if (numel(files_E)~=numel(files_H))
            error 'number of E/H planes mismatch!'
        end
        for fn = 1:numel(files_E)
            nf2ff_xml.Planes{end+1}.ATTRIBUTE.E_Field = files_E(fn).name;
            nf2ff_xml.Planes{end}.ATTRIBUTE.H_Field = files_H(fn).name;
        end
    end
end

nf2ff_xml.ATTRIBUTE.freq = freq;
nf2ff_xml.theta = theta;
nf2ff_xml.phi = phi;

nf2ff.xml = [Sim_Path '' filesep '' filename '.xml'];
nf2ff.hdf5 = [Sim_Path '' filesep '' nf2ff_xml.ATTRIBUTE.Outfile];

% create nf2ff structure
struct_2_xml(nf2ff.xml,nf2ff_xml,'nf2ff');

m_filename = mfilename('fullpath');
dir_name = fileparts( m_filename );

if isunix
    nf2ff_bin = searchBinary('nf2ff', ...
    {[dir_name filesep '..' filesep 'nf2ff' filesep], ...
     [dir_name filesep '..' filesep '..' filesep '..' filesep 'bin' filesep]}, 0);
else
    nf2ff_bin = searchBinary('nf2ff.exe',[dir_name filesep '..' filesep], 0);
end

if ((exist(nf2ff.hdf5,'file') && (mode==0)) || (mode==2))
    disp('CalcNF2FF: Reading nf2ff data only...')
    nf2ff = ReadNF2FF(nf2ff);

    % verify read data
    if ( (vectorEqual(nf2ff.freq,freq)==0) || (vectorEqual(nf2ff.theta,theta)==0) || (vectorEqual(nf2ff.phi,phi)==0) )
        error('openEMS:CalcNF2FF','data mismatch between read and requested data --> recalculate nf2ff --> Set Mode to 1 ');
    end
    return;
end

savePath = pwd;
cd(Sim_Path);

try
    if (isempty(nf2ff_bin))
        error('openEMS:CalcNF2FF','nf2ff binary not found!');
    end
    if isunix
        % remove LD_LIBRARY_PATH set by matlab
        system(['export LD_LIBRARY_PATH=; ' nf2ff_bin ' ' filename '.xml']);
    else
        system([nf2ff_bin ' ' filename '.xml']);
    end
    nf2ff.hdf5;
    cd(savePath);
catch
    cd(savePath);
    error 'CalcNF2FF: failed'
end

nf2ff = ReadNF2FF(nf2ff);

% verify read data
if ( (vectorEqual(nf2ff.freq,freq)==0) || (vectorEqual(nf2ff.theta,theta)==0) || (vectorEqual(nf2ff.phi,phi)==0) )
    error('openEMS:CalcNF2FF','data mismatch between read and requested data --> THIS SHOULD NOT HAPPEN!');
end

function equal = vectorEqual(v1, v2, acc)
if (nargin<3)
    acc = 1e-6;
end

equal = 0;
if numel(v1)~=numel(v2)
    return;
end

if sum(abs((v1(:)-v2(:))/v1(:)) > acc)>0
    return;
end
equal = 1;
return
