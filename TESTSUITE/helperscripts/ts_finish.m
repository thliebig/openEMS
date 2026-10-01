function pass = ts_finish(opt, name, pass, Sim_Path)
%pass = ts_finish(opt, name, pass, Sim_Path)
%
% Common end of a test: print the overall verdict, drop the simulation folder
% of a passed test if opt.Cleanup is set and raise an error on failure if
% opt.StopIfFailed is set.
%
% Sim_Path is optional; tests that do not simulate leave it out.
%
% openEMS testsuite
% -----------------
%
% See also ts_check, ts_check_rel, ts_options

pass = ~isempty(pass) && all(pass(:));

if pass
    fprintf('  ==> %s: PASS\n', name);
else
    fprintf('  ==> %s: FAILED\n', name);
end
if isOctave
    fflush(stdout);
end

if pass && opt.Cleanup && (nargin > 3) && ~isempty(Sim_Path) && exist(Sim_Path,'dir')
    rmdir(Sim_Path, 's');
end

if ~pass && opt.StopIfFailed
    error('openEMS:TESTSUITE', '%s failed', name);
end
