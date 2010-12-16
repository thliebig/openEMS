function UI = ReadUI(files, path, freq)
% function UI = ReadUI(files, path, freq)
%
% read current and voltages from multiple files found in path
%
% returns voltages/currents in time and frequency-domain
% 
% remarks on the frequency-domain:
% - all signals are assumed to start at t=0
% - currents that e.g. start at t = +delta_t/2 will be phase shifted by
%    exp(-j*w*t(1))
%
% parameter:
% freq (optional):
%     frequency-domain values will be calculated according to 'freq'
%     if 'freq' is not given, a FFT will be used
%
% e.g.
% U = ReadUI({'ut1_1','ut1_2'},'tmp' );
% I = ReadUI('it1'            ,'tmp',[0.5e9 1e9 1.5e9]);
% 
% openEMS matlab interface
% -----------------------
% author: Thorsten Liebig

if (nargin<2)
    path ='';
end

if (ischar(files))
    filenames{1}=files;
else
    filenames=files;
end

UI.TD = {};
UI.FD = {};
for n=1:numel(filenames)
    tmp = load( fullfile(path,filenames{n}) );
    t = tmp(:,1)';
    val = tmp(:,2)';
    
    UI.TD{n}.t = t;
    UI.TD{n}.val = val;

    if (numel(tmp(1,:))>2)
        UI.TD{n}.additional = tmp(:,3:end)';
    end    
    
    if (nargin<3) || isempty(freq)
        [UI.FD{n}.f,UI.FD{n}.val] = FFT_time2freq( t,val );
    else
        UI.FD{n}.f = freq;
        UI.FD{n}.val = DFT_time2freq( t, val, freq );
    end
    
    %correct phase error for time-shifted signals
    UI.FD{n}.val = UI.FD{n}.val .* exp(-1j*2*pi*UI.FD{n}.f * t(1)); 
end
