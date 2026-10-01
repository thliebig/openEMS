%run_testsuite  run the openEMS Octave/Matlab test suite
%
% Every test is a function in one of the group folders next to this script
% and returns a pass/fail verdict:
%
%     pass = <test>(<key>, <value>, ...)        see ts_options
%
% This script runs them all, reports PASS/FAIL/ERROR per test and prints a
% summary. A failing test neither aborts the run nor hides the tests behind
% it.
%
% usage:
%
%     octave --no-gui run_testsuite.m [<options>] [<name>...]
%
% options:
%
%   --engine=<e>    run with this engine, e.g. --engine=basic. By default
%                   openEMS picks the fastest engine it has; comparing the
%                   engines is the job of enginetests/, not of every test.
%   --all-engines   run every simulating test with every engine (slow)
%   --keep          keep the simulation folder of passed tests too
%   --list          list the tests without running them
%   <name>          run only the tests whose group/name contains <name>
%
% Inside Octave or Matlab just call "run_testsuite" for a plain default run.
% A script cannot take arguments, so to pass the options above from a session,
% set them first; they apply to the next run only:
%
%     run_testsuite_args = {'--list'};
%     run_testsuite
%
% Or call a single test directly, with plots:
%
%     cd combinedtests; Coax('Plots',1)
%
% openEMS testsuite
% -----------------
%
% See also ts_options, ts_check, ts_finish

% ---------------------------------------------------------------- environment
ts_folder = fileparts(mfilename('fullpath'));
addpath(fullfile(ts_folder, 'helperscripts'));

if ~exist('InitFDTD','file') || ~exist('InitCSX','file') || ~exist('isOctave','file')
    error('openEMS:TESTSUITE', ...
        ['the openEMS and CSXCAD matlab interfaces are not on the path:\n' ...
         '  addpath(''<prefix>/share/openEMS/matlab'')\n' ...
         '  addpath(''<prefix>/share/CSXCAD/matlab'')']);
end

if isOctave()
    confirm_recursive_rmdir(0);
    page_screen_output(0);      % do not buffer output
    page_output_immediately(1); % do not buffer output
end

% -------------------------------------------------------------------- options
ts_engines = {''};      % '' -> let openEMS choose the fastest engine
ts_keep    = 0;
ts_list    = 0;
ts_filter  = {};

% groups that do their own engine handling and are never swept over engines
ts_engine_independent = {'unittests', 'enginetests'};

% all known engines, in the order the engine sweep uses them
ts_all_engines = {'--engine=basic', '--engine=sse', ...
                  '--engine=sse-compressed', '--engine=multithreaded'};

% Only "octave run_testsuite.m [options]" puts this script's options into
% argv(); in a session argv() holds Octave's own command line instead (--gui and
% friends), so it must not be parsed there.
ts_from_cmdline = isOctave() && ...
                  ~isempty(regexp(program_invocation_name, '\.m$', 'once'));

ts_args = {};
if ts_from_cmdline
    ts_args = argv();
elseif exist('run_testsuite_args', 'var')
    ts_args = run_testsuite_args;    % the way to pass options from a session
    clear run_testsuite_args         % one run, so it cannot linger unnoticed
    if ischar(ts_args)
        ts_args = {ts_args};
    end
end
for ts_n = 1:numel(ts_args)
    ts_arg = ts_args{ts_n};
    if strncmp(ts_arg, '--engine=', 9)
        % openEMS silently falls back to its default on an unknown engine
        % name, so a typo here would quietly test something else
        if ~any(strcmp(ts_arg, ts_all_engines)) && ~strcmp(ts_arg, '--engine=fastest')
            error('openEMS:TESTSUITE', 'unknown engine: %s\n  known: %s', ...
                  ts_arg(10:end), strjoin(strrep(ts_all_engines,'--engine=',''), ', '));
        end
        ts_engines = {ts_arg};
    elseif strcmp(ts_arg, '--all-engines')
        ts_engines = ts_all_engines;
    elseif strcmp(ts_arg, '--keep')
        ts_keep = 1;
    elseif strcmp(ts_arg, '--list')
        ts_list = 1;
    elseif strncmp(ts_arg, '-', 1)
        error('openEMS:TESTSUITE', 'unknown option: %s', ts_arg);
    else
        ts_filter{end+1} = ts_arg;
    end
end

% ---------------------------------------------------------- collect the tests
% cheap groups first, so a broken installation shows up in seconds
ts_group_order = {'unittests', 'probes', 'combinedtests', 'enginetests'};

ts_entries = dir(ts_folder);
ts_found = {};
for ts_n = 1:numel(ts_entries)
    if ~ts_entries(ts_n).isdir, continue; end
    ts_name = ts_entries(ts_n).name;
    if ts_name(1)=='.' || strcmp(ts_name,'helperscripts') || strncmp(ts_name,'tmp',3)
        continue
    end
    ts_found{end+1} = ts_name;
end
% known groups in the order above, anything new appended alphabetically
ts_groups = {};
for ts_n = 1:numel(ts_group_order)
    if any(strcmp(ts_group_order{ts_n}, ts_found))
        ts_groups{end+1} = ts_group_order{ts_n};
    end
end
ts_groups = [ts_groups, sort(setdiff(ts_found, ts_group_order))];

ts_tests = {};   % {group, name} per test
for ts_g = 1:numel(ts_groups)
    ts_scripts = dir(fullfile(ts_folder, ts_groups{ts_g}, '*.m'));
    for ts_s = 1:numel(ts_scripts)
        if ts_scripts(ts_s).isdir, continue; end
        [~, ts_name] = fileparts(ts_scripts(ts_s).name);
        ts_id = [ts_groups{ts_g} '/' ts_name];
        if ~isempty(ts_filter) && ...
           ~any(cellfun(@(p) ~isempty(strfind(ts_id, p)), ts_filter))
            continue
        end
        ts_tests{end+1} = {ts_groups{ts_g}, ts_name};
    end
end

if isempty(ts_tests)
    error('openEMS:TESTSUITE', 'no test matches the given name(s)');
end

if ts_list
    fprintf('*** openEMS testsuite -- available tests:\n');
    for ts_n = 1:numel(ts_tests)
        fprintf('  %s/%s\n', ts_tests{ts_n}{1}, ts_tests{ts_n}{2});
    end
    return
end

% ------------------------------------------------------------------- run them
ts_results = {};   % {id, engine, status, time, message} per run
ts_total   = tic();

for ts_e = 1:numel(ts_engines)
    ts_engine = ts_engines{ts_e};
    if isempty(ts_engine)
        ts_engine_label = 'default (fastest available)';
    else
        ts_engine_label = ts_engine;
    end
    fprintf('\n*** %s openEMS testsuite started -- engine: %s\n\n', ...
            datestr(now), ts_engine_label);

    for ts_n = 1:numel(ts_tests)
        ts_group = ts_tests{ts_n}{1};
        ts_name  = ts_tests{ts_n}{2};
        ts_id    = [ts_group '/' ts_name];

        if ts_e > 1 && any(strcmp(ts_group, ts_engine_independent))
            continue   % already covered by the first engine
        end

        fprintf('[%2d/%2d] %s\n', ts_n, numel(ts_tests), ts_id);
        if isOctave(), fflush(stdout); end

        ts_oldpwd = pwd();
        ts_tic    = tic();
        ts_msg    = '';
        try
            cd(fullfile(ts_folder, ts_group));
            ts_pass = feval(ts_name, 'openEMS_opts', ts_engine, 'Plots', 0, ...
                            'Silent', 1, 'Verbose', 0, 'Cleanup', ~ts_keep, ...
                            'StopIfFailed', 0);
            if ts_pass
                ts_status = 'PASS';
            else
                ts_status = 'FAIL';
            end
        catch ts_err
            ts_status = 'ERROR';
            ts_msg    = ts_err.message;
            fprintf('  ==> %s: ERROR -- %s\n', ts_name, ts_msg);
        end
        cd(ts_oldpwd);

        ts_results{end+1} = {ts_id, ts_engine_label, ts_status, toc(ts_tic), ts_msg};
    end
end

% -------------------------------------------------------------------- summary
ts_n_pass  = 0;
ts_n_fail  = 0;
ts_n_error = 0;

fprintf('\n*** openEMS testsuite summary\n\n');
for ts_n = 1:numel(ts_results)
    ts_r = ts_results{ts_n};
    fprintf('  %-7s %-32s %7.1f s', ts_r{3}, ts_r{1}, ts_r{4});
    if numel(ts_engines) > 1
        fprintf('  [%s]', ts_r{2});
    end
    fprintf('\n');
    if strcmp(ts_r{3},'PASS')
        ts_n_pass = ts_n_pass + 1;
    elseif strcmp(ts_r{3},'FAIL')
        ts_n_fail = ts_n_fail + 1;
    else
        ts_n_error = ts_n_error + 1;
    end
end

fprintf('\n  %d passed, %d failed, %d errors in %.1f s\n', ...
        ts_n_pass, ts_n_fail, ts_n_error, toc(ts_total));

if ts_n_fail + ts_n_error > 0
    fprintf('\n*** TESTSUITE FAILED\n');
else
    fprintf('\n*** ALL TESTS PASSED\n');
end
if isOctave(), fflush(stdout); end

% started as "octave run_testsuite.m": report the verdict as an exit status
if ts_from_cmdline
    if ts_n_fail + ts_n_error > 0
        exit(1);
    else
        exit(0);
    end
end
