function Sim_Path = ts_sim_path(test_fullpath)
%Sim_Path = ts_sim_path(test_fullpath)
%
% Return an empty simulation folder for a test, created next to the test itself
% as tmp_<testname>. Call it as
%
%     Sim_Path = ts_sim_path(mfilename('fullpath'));
%
% so the test works no matter what the current directory is, and so a run never
% sees the output of the previous one. Failing to clear the folder is an error
% rather than something a test could silently pass on stale results with.
%
% openEMS testsuite
% -----------------
%
% See also ts_finish, ts_options

[folder, name] = fileparts(test_fullpath);
Sim_Path = fullfile(folder, ['tmp_' name]);

if exist(Sim_Path, 'dir')
    if isOctave()
        % do not rely on the caller having turned the prompt off: a test run
        % from the Octave prompt has to clear its folder just as reliably
        old_confirm = confirm_recursive_rmdir(false);
        [ok, msg] = rmdir(Sim_Path, 's');
        confirm_recursive_rmdir(old_confirm);
    else
        [ok, msg] = rmdir(Sim_Path, 's');
    end
    if ~ok
        error('openEMS:TESTSUITE', 'cannot clear %s: %s', Sim_Path, msg);
    end
end

[ok, msg] = mkdir(Sim_Path);
if ~ok
    error('openEMS:TESTSUITE', 'cannot create %s: %s', Sim_Path, msg);
end
