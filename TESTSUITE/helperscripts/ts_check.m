function pass = ts_check(label, pass, varargin)
%pass = ts_check(label, pass, <fmt>, <args>...)
%
% Report a single check of a test and pass its verdict through, so checks can
% be collected with &&:
%
%     ok = ts_check('char. impedance', abs(Z-Z0)<tol, '%.2f Ohm', Z);
%
% The optional format string and arguments are printed behind the label and
% should carry the measured value and the tolerance -- that is what one needs
% to judge a failure from a log file.
%
% openEMS testsuite
% -----------------
%
% See also ts_check_rel, ts_finish, ts_options

pass = ~isempty(pass) && all(pass(:));

if nargin > 2
    detail = [' -- ' sprintf(varargin{:})];
else
    detail = '';
end

if pass
    tag = 'PASS';
else
    tag = 'FAIL';
end

fprintf('    [%s] %s%s\n', tag, label, detail);
if isOctave
    fflush(stdout);
end
