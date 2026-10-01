function [pass, err] = ts_check_rel(label, value, reference, tol)
%[pass, err] = ts_check_rel(label, value, reference, tol)
%
% Check every element of value against reference within the relative
% tolerance tol and report the result. reference is either a scalar or has as
% many elements as value. Returns the largest relative deviation found.
%
% openEMS testsuite
% -----------------
%
% See also ts_check, ts_finish, ts_options

value     = value(:);
reference = reference(:);

if numel(reference) == 1
    reference = repmat(reference, size(value));
elseif numel(reference) ~= numel(value)
    error('openEMS:TESTSUITE','value and reference must have the same size');
end

if any(reference == 0)
    error('openEMS:TESTSUITE','relative check against a zero reference');
end

err  = max(abs((value - reference) ./ reference));
pass = ts_check(label, err <= tol, ...
                'max. rel. deviation %.3g (tolerance %.3g)', err, tol);
