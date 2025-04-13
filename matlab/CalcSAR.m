function CalcSAR(sar_fn, sar_out, varargin)
% function CalcSAR(sar_fn, sar_out, varargin)
%
% Calculate the SAR (specific absorption rate)
%
% Note: No averaging method like IEEE_62704 are (yet) validated
%       according to the IEC/IEEE-62704-1!
%
% parameter:
% - sar_fn:   hdf5 file with SAR raw data as created by a dumbox
% - sar_out:  hdf5 output file with SAR results
%
% optional parameter:
% - 'mass':   averaging mass in g (default is no averaging)
% - 'method': IEEE_C95_3, IEEE_62704 or SIMPLE (default)
% - 'verbose': verbose output
%
% See also: AddDump
%
% openEMS matlab interface
% -----------------------
% author: Thorsten Liebig, 2025

bin_args = ' ';

for n=1:2:numel(varargin)
    if (strcmp(varargin{n},'verbose'))
        bin_args = [bin_args ' -v'];
    elseif (strcmp(varargin{n},'mass'))
        bin_args = [bin_args ' --mass ' num2str(varargin{n+1})];
    elseif (strcmp(varargin{n},'method'))
        bin_args = [bin_args ' --method ' varargin{n+1}];
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
        error('openEMS:CalcSAR','nf2ff binary not found!');
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

