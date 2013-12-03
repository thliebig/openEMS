function [CSX nf2ff] = CreateNF2FFBox(CSX, name, start, stop, varargin)
% function [CSX nf2ff] = CreateNF2FFBox(CSX, name, start, stop, varargin)
%
% create the dump boxes needed for the near field to far field transformation
% 
% input:
%   name:       name of this nf2ff box
%   start/stop: start/stop coordinates for the nf2ff box (this box has to
%               enclose all radiating structures!)
% optional inputs:
%   'Directions': enable/disable specific directions, e.g.
%                 'Directions',[1 1 0 0 1 1]
%                   -> disable nf2ff in +/-y direction
%   'Frequency': dump nf2ff in frequency domain, this will save disk-space
%                but is less flexible, since only this frequencies can be
%                used for the nf2ff calculations by CalcNF2FF
%                See also AddDump for more information
%   'OptResolution': specify a dump resolution, this will save disk-space
%                    See also AddDump for more information
%                    e.g.: 'OptResolution', c0/max_freq/unit/15
%
% example:
%   see Tutorials/Simple_Patch_Antenna.m
%   see Tutorials/Helical_Antenna.m
% 
% See also CalcNF2FF
%
% (C) 2010 Sebastian Held <sebastian.held@gmx.de>
% (C) 2010-2012 Thorsten Liebig <thorsten.liebig@gmx.de>

if (nargin<5)
    directions = ones(6,1);
end

directions = ones(6,1);
add_args = {};
dump_type = 0;
dump_mode = 1;

for n=1:numel(varargin)/2
    if strcmp(varargin{2*n-1},'Frequency')
        add_args = {add_args{:}, 'Frequency', varargin{2*n}};
        dump_type = 10;
    elseif strcmp(varargin{2*n-1},'Directions')
        directions=varargin{2*n};
    else
        add_args = {add_args{:}, varargin{2*n-1}, varargin{2*n}};
    end
end

nf2ff.name = name;
nf2ff.filenames_E = {[name '_E_xn'],[name '_E_xp'],[name '_E_yn'],[name '_E_yp'],[name '_E_zn'],[name '_E_zp']};
nf2ff.filenames_H = {[name '_H_xn'],[name '_H_xp'],[name '_H_yn'],[name '_H_yp'],[name '_H_zn'],[name '_H_zp']};
nf2ff.directions = directions;

if (isfield(CSX,'ATTRIBUTE'))
    if (isfield(CSX.ATTRIBUTE,'CoordSystem'))
        nf2ff.CoordSystem = CSX.ATTRIBUTE.CoordSystem;
    end
    if (isfield(CSX,'BackgroundMaterial'))
        if (isfield(CSX.ATTRIBUTE,'Epsilon'))
            nf2ff.Eps_r = CSX.ATTRIBUTE.BG_epsR;
        end
        if (isfield(CSX.ATTRIBUTE,'Mue'))
            nf2ff.Mue_r = CSX.ATTRIBUTE.BG_mueR;
        end
    end
end

for nd = 1:3
    pos = 2*nd-1;
    if (directions(pos))
        l_start = start;
        l_stop = stop;
        l_stop(nd) = start(nd);
        CSX = AddBox( AddDump(CSX,nf2ff.filenames_E{pos},'DumpType',dump_type,'DumpMode',dump_mode,'FileType',1,add_args{:}), nf2ff.filenames_E{pos}, 0, l_start, l_stop );
        CSX = AddBox( AddDump(CSX,nf2ff.filenames_H{pos},'DumpType',dump_type+1,'DumpMode',dump_mode,'FileType',1,add_args{:}), nf2ff.filenames_H{pos}, 0, l_start, l_stop );
    else
        nf2ff.filenames_E{pos}='';
        nf2ff.filenames_H{pos}='';
    end
    pos = 2*nd;
    if (directions(pos))
        l_start = start;
        l_stop = stop;
        l_start(nd) = stop(nd);
        CSX = AddBox( AddDump(CSX,nf2ff.filenames_E{pos},'DumpType',dump_type,'DumpMode',dump_mode,'FileType',1,add_args{:}), nf2ff.filenames_E{pos}, 0, l_start, l_stop );
        CSX = AddBox( AddDump(CSX,nf2ff.filenames_H{pos},'DumpType',dump_type+1,'DumpMode',dump_mode,'FileType',1,add_args{:}), nf2ff.filenames_H{pos}, 0, l_start, l_stop );
    else
        nf2ff.filenames_E{pos}='';
        nf2ff.filenames_H{pos}='';
    end
end

