% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2026 openEMS contributors
function snapshot_path = RequestDFTSnapshot(sim_path, request_id)
% Request an incomplete E/H DFT export from a running --dft-snapshots solver.
% Call from another session. Use one producer and one outstanding request.
% Wait for snapshot_path/snapshot.xml before reading any of its HDF5 files.
% A sibling <request_id>.error reports failure. This does not stop the solver
% and is not a restart checkpoint. The run must use sim_path as its cwd.
if isempty(regexp(request_id, '^[A-Za-z0-9_-]{1,64}$', 'once'))
    error('Invalid snapshot ID: use 1-64 ASCII letters/digits/-/_');
end
snapshot_path = fullfile(sim_path, 'dft-snapshots', request_id);
if exist(snapshot_path, 'dir') || exist([snapshot_path '.error'], 'file')
    error('Snapshot ID already used');
end
request = fullfile(sim_path, 'DFT_SNAPSHOT');
if exist(request, 'file') || exist([request '.active'], 'file')
    error('A snapshot request is already pending');
end
temporary = tempname(sim_path);
fid = fopen(temporary, 'w');
if fid < 0
    error('Cannot create snapshot request');
end
fprintf(fid, '%s\n', request_id);
fclose(fid);
[ok, message] = movefile(temporary, request);
if ~ok
    delete(temporary);
    error('Cannot publish snapshot request: %s', message);
end
end
