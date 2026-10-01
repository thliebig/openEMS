function opt = ts_options(varargin)
%opt = ts_options(<key>, <value>, ...)
%
% Parse the options every test of this suite accepts and fill in the
% defaults. A test starts with
%
%     opt = ts_options(varargin{:});
%
% and then honours opt.Plots, opt.Silent and so on.
%
% - 'Engine':       additional openEMS command line options, usually the
%                   engine selection, e.g. '--engine=basic'. The default ''
%                   lets openEMS pick the fastest available engine.
% - 'Plots':        0/1, draw the diagnostic plots of the test. The default is
%                   on, unless there is no graphics toolkit to draw with -- a
%                   single test run over ssh should report its result, not die
%                   in figure().
% - 'Silent':       0/1, suppress the openEMS console output (default 0)
% - 'Verbose':      0/1, print additional details while running (default 1)
% - 'Cleanup':      0/1, remove the simulation folder of a passed test
%                   (default 1)
% - 'StopIfFailed': 0/1, raise an error if the test failed (default 1)
%
% run_testsuite calls every test with 'Plots',0, 'Silent',1, 'Verbose',0 and
% 'StopIfFailed',0, so that one failing test neither blocks nor aborts the
% run.
%
% openEMS testsuite
% -----------------
%
% See also ts_check, ts_check_rel, ts_finish, run_testsuite

opt.openEMS_opts     = '';
opt.ExactEndCriteria = 1;
opt.Plots            = can_plot();
opt.Silent           = 0;
opt.Verbose          = 1;
opt.Cleanup          = 1;
opt.StopIfFailed     = 1;

if mod(numel(varargin),2) ~= 0
    error('openEMS:TESTSUITE','test options must be given as key/value pairs');
end

known = fieldnames(opt);
for n=1:2:numel(varargin)
    key = varargin{n};
    if ~ischar(key)
        error('openEMS:TESTSUITE','test option %d is not a key', (n+1)/2);
    end
    idx = find(strcmpi(known, key));
    if isempty(idx)
        error('openEMS:TESTSUITE','unknown test option: %s', key);
    end
    opt.(known{idx}) = varargin{n+1};
end

% a test suite wants a stopping point that does not depend on the machine
if opt.ExactEndCriteria
    opt.openEMS_opts = strtrim([opt.openEMS_opts ' --exact-endcriteria']);
end

end


function tf = can_plot()
tf = true;
if isOctave()
    tf = ~isempty(available_graphics_toolkits());
end
end
