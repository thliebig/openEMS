function CalcSAR(sar_fn, sar_out, varargin)
% function CalcSAR(sar_fn, sar_out, varargin)
%
% Calculate the SAR (specific absorption rate)
%
% Note: No averaging method like IEEE_62704 are (yet) validated
%       according to the IEC/IEEE-62704-1!
%
% parameter:
% - sar_fn:   hdf5 file with SAR raw data as created by a dump box
% - sar_out:  hdf5 output file with SAR results
%
% optional parameter:
% - 'mass':       averaging mass in g (default is no averaging)
% - 'method':     IEEE_C95_3, IEEE_62704 or SIMPLE (default)
%                 All three are real cubical mass averaging methods and
%                 differ only in how strictly an averaging cube has to be
%                 valid; SIMPLE accepts any cube reaching the target mass,
%                 the IEEE methods reject cubes at a boundary and fill
%                 those cells in from a neighbour. Has no effect if 'mass'
%                 is 0, i.e. local SAR without any averaging.
% - 'autoRange':  limit the calculation to the cells within N dB of the peak
%                 local SAR, plus a padding of about one averaging cube. The
%                 result is returned on this reduced mesh. This is a speedup
%                 and not a guarantee to find the global peak, do not use it
%                 for standard compliance work.
% - 'numThreads': number of worker threads (default: all available)
% - 'progress':   show progress output
% - 'verbose':    verbose output
%
% See also: AddDump
%
% openEMS matlab interface
% -----------------------
% author: Thorsten Liebig, 2025

bin_args = ' --legacyHDF5Dumps';

for n=1:2:numel(varargin)
    if (strcmp(varargin{n},'verbose'))
        bin_args = [bin_args ' -v'];
    elseif (strcmp(varargin{n},'progress'))
        bin_args = [bin_args ' -p'];
    elseif (strcmp(varargin{n},'mass'))
        bin_args = [bin_args ' --mass ' num2str(varargin{n+1})];
    elseif (strcmp(varargin{n},'method'))
        bin_args = [bin_args ' --method ' varargin{n+1}];
    elseif (strcmp(varargin{n},'autoRange'))
        bin_args = [bin_args ' --autorange ' num2str(varargin{n+1})];
    elseif (strcmp(varargin{n},'numThreads'))
        bin_args = [bin_args ' --numThreads ' num2str(varargin{n+1})];
    else
        error(['Invalid argument: ' varargin{n}]);
    end
end


m_filename = mfilename('fullpath');
dir_name = fileparts( m_filename );

if isunix
    sar_calc_bin = searchBinary('sar_calc', ...
    {[dir_name filesep '..' filesep '..' filesep '..' filesep 'bin' filesep]}, 0);
else
    sar_calc_bin = searchBinary('sar_calc.exe', [dir_name filesep '..' filesep], 0);
end

try
    if (isempty(sar_calc_bin))
        error('openEMS:CalcSAR','sar_calc binary not found!');
    end
    cmd = [sar_calc_bin ' -i ' sar_fn ' -o ' sar_out bin_args];
    if isunix
        % remove LD_LIBRARY_PATH set by matlab
        system(['export LD_LIBRARY_PATH=; ' cmd]);
    else
        system(cmd);
    end
catch
    error 'CalcSAR: failed'
end

