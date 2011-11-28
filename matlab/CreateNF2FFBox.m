function [CSX nf2ff] = CreateNF2FFBox(CSX, name, start, stop, directions)
% function [CSX nf2ff] = CreateNF2FFBox(CSX, name, start, stop, directions)
%
% create the dump boxes needed for the near field to far field transformation
% 
% input:
%   name:       name of this nf2ff box
%   start/stop: start/stop coordinates for the nf2ff box (this box has to
%               enclose all radiating structures!)
%
% example:
%   see examples/NF2FF/infDipol.m
% 
% See also AnalyzeNF2FF
%
% (C) 2010 Sebastian Held <sebastian.held@gmx.de>
% (C) 2010, 2011 Thorsten Liebig <thorsten.liebig@gmx.de>

if (nargin<5)
    directions = ones(6,1);
end

nf2ff.filenames_E = {[name '_Et_xn'],[name '_Et_xp'],[name '_Et_yn'],[name '_Et_yp'],[name '_Et_zn'],[name '_Et_zp']};
nf2ff.filenames_H = {[name '_Ht_xn'],[name '_Ht_xp'],[name '_Ht_yn'],[name '_Ht_yp'],[name '_Ht_zn'],[name '_Ht_zp']};
nf2ff.directions = directions;

if (directions(1))
    CSX = AddBox( AddDump(CSX,nf2ff.filenames_E{1},'DumpType',0,'DumpMode',2,'FileType',1), nf2ff.filenames_E{1}, 0, start, [start(1) stop(2) stop(3)] );
    CSX = AddBox( AddDump(CSX,nf2ff.filenames_H{1},'DumpType',1,'DumpMode',2,'FileType',1), nf2ff.filenames_H{1}, 0, start, [start(1) stop(2) stop(3)] );
end

if (directions(2))
    CSX = AddBox( AddDump(CSX,nf2ff.filenames_E{2},'DumpType',0,'DumpMode',2,'FileType',1), nf2ff.filenames_E{2}, 0, [stop(1) start(2) start(3)], stop );
    CSX = AddBox( AddDump(CSX,nf2ff.filenames_H{2},'DumpType',1,'DumpMode',2,'FileType',1), nf2ff.filenames_H{2}, 0, [stop(1) start(2) start(3)], stop );
end

if (directions(3))
    CSX = AddBox( AddDump(CSX,nf2ff.filenames_E{3},'DumpType',0,'DumpMode',2,'FileType',1), nf2ff.filenames_E{3}, 0, start, [stop(1) start(2) stop(3)] );
    CSX = AddBox( AddDump(CSX,nf2ff.filenames_H{3},'DumpType',1,'DumpMode',2,'FileType',1), nf2ff.filenames_H{3}, 0, start, [stop(1) start(2) stop(3)] );
end

if (directions(4))
    CSX = AddBox( AddDump(CSX,nf2ff.filenames_E{4},'DumpType',0,'DumpMode',2,'FileType',1), nf2ff.filenames_E{4}, 0, [start(1) stop(2) start(3)], stop );
    CSX = AddBox( AddDump(CSX,nf2ff.filenames_H{4},'DumpType',1,'DumpMode',2,'FileType',1), nf2ff.filenames_H{4}, 0, [start(1) stop(2) start(3)], stop );
end

if (directions(5))
    CSX = AddBox( AddDump(CSX,nf2ff.filenames_E{5},'DumpType',0,'DumpMode',2,'FileType',1), nf2ff.filenames_E{5}, 0, start, [stop(1) stop(2) start(3)] );
    CSX = AddBox( AddDump(CSX,nf2ff.filenames_H{5},'DumpType',1,'DumpMode',2,'FileType',1), nf2ff.filenames_H{5}, 0, start, [stop(1) stop(2) start(3)] );
end

if (directions(6))
    CSX = AddBox( AddDump(CSX,nf2ff.filenames_E{6},'DumpType',0,'DumpMode',2,'FileType',1), nf2ff.filenames_E{6}, 0, [start(1) start(2) stop(3)], stop );
    CSX = AddBox( AddDump(CSX,nf2ff.filenames_H{6},'DumpType',1,'DumpMode',2,'FileType',1), nf2ff.filenames_H{6}, 0, [start(1) start(2) stop(3)], stop );
end
