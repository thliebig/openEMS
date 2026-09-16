function ensure_oct_file(fname)
% ensure_oct_file(fname)
%
% internal function for openEMS to make sure the oct-file <fname> can be
% used, running "setup" if it cannot.
%
% An oct-file belongs to the Octave version it was compiled with. After an
% Octave update the old file is still found by exist(), but loading it fails
% with "Incompatible version or missing dependency?", so looking for the file
% alone is not enough to decide whether "setup" has to run.
%
% See also: setup
%
% openEMS Matlab/Octave interface
% -----------------------
% author: Thorsten Liebig, 2026

if (oct_file_usable(fname))
    return
end

warning('openEMS:ensure_oct_file', ['oct-file "' fname '" is missing or was built for another Octave version, trying to run "setup"']);
try
    setup
catch err
    error('openEMS:ensure_oct_file', ['running "setup" failed: ' err.message]);
end

% drop the failed load, so that the freshly compiled file is picked up
clear('-f', fname);
if (~oct_file_usable(fname))
    error('openEMS:ensure_oct_file', ['oct-file "' fname '" is still not usable after running "setup"']);
end

end

function usable = oct_file_usable(fname)
% Calling the function without arguments loads the oct-file and is answered
% with its usage message; a file built for another Octave version fails to
% load instead. There is no way to load it without calling it.
usable = false;
if (exist(fname) == 0)
    return
end
usable = true;
try
    feval(fname);
catch err
    usable = isempty(strfind(err.message, 'failed to load'));
end

end
