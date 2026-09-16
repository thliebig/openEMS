function h5writemode(filename, x, y, Vx, Vy)
% h5writemode(filename, x, y, Vx, Vy)
%
% Write a 2-D mode field profile to an HDF5 file for use with
% AddWaveGuidePort ('E_WG_file' / 'H_WG_file').
%
% Parameters
% ----------
%   filename : output file path (overwritten if it exists)
%   x, y     : coordinate vectors of length nx and ny
%   Vx, Vy   : field component matrices (nx-by-ny), where
%              Vx(i,j) is the field value at coordinate (x(i), y(j))
%
% HDF5 layout written:
%   /x   1-D dataset, length nx
%   /y   1-D dataset, length ny
%   /Vx  2-D dataset, shape nx x ny
%   /Vy  2-D dataset, shape nx x ny
%   root attribute Version = 1.0
%
% On Octave the oct-file h5writemode_octave (compiled by setup.m) is used.
% On Matlab  the built-in h5create / h5write / h5writeatt are used.
%
% See also AddWaveGuidePort, setup

if isOctave()
    ensure_oct_file('h5writemode_octave');
    % Delegate to the oct-file compiled by setup.m.
    % Ensure x and y are row vectors so matrix_value() gives (1 x n).
    h5writemode_octave(filename, x(:)', y(:)', Vx, Vy);
else
    % Matlab: h5write reverses array dimensions (column-major -> row-major).
    % A Matlab (ny x nx) matrix becomes an HDF5 (nx x ny) dataset, which is
    % exactly what CSModeData expects. Write the transposed Vx/Vy.
    nx = numel(x);
    ny = numel(y);
    if exist(filename, 'file'), delete(filename); end
    h5create(filename, '/x',  nx);
    h5write( filename, '/x',  x(:)');
    h5create(filename, '/y',  ny);
    h5write( filename, '/y',  y(:)');
    h5create(filename, '/Vx', [ny, nx]);
    h5write( filename, '/Vx', Vx.');
    h5create(filename, '/Vy', [ny, nx]);
    h5write( filename, '/Vy', Vy.');
    h5writeatt(filename, '/', 'Version', 1.0);
end
