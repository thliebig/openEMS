function UI = ReadUI(files, path, freq, varargin)
% function UI = ReadUI(files, path, freq, varargin)
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
% optional parameter:
% freq: frequency-domain values will be calculated according to 'freq'
%       if 'freq' is not given, a (zero padded) FFT will be used
%
% optional key,value pairs:
% 'AR'  : auto-regressive model to improve FD accuracy
%         values: order to use within an AR model or 'auto'
%
% % examples:
% U = ReadUI({'ut1_1','ut1_2'},'tmp' );
% I = ReadUI('it1'            ,'tmp',[0.5e9 1e9 1.5e9]);
% 
% % using the auto-regressive model
% U = ReadUI('port_ut1' , 'tmp', 'AR', 'auto');
% 
% openEMS matlab interface
% -----------------------
% author: Thorsten Liebig
%
% See also DFT_time2freq, AR_estimate

if (nargin<2)
    path ='';
end

AR_order = 0;
SignalType = 'pulse';

for n=1:2:numel(varargin)
    if (strcmp(varargin{n},'AR')==1)
        AR_order =  varargin{n+1};
    elseif strcmpi(varargin{n},'SignalType')
        SignalType = varargin{n+1};
    else
        warning('CSXCAD:ReadUI', ['"' varargin{n} '" is an unknown argument']);
    end
end

if strcmpi(SignalType,'periodic') && AR_order>0
    error 'auto-regressive model not compatible with periodic signals'
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
        if strcmpi(SignalType,'periodic')
           warning 'ReadUI: periodic signal type not supported by FFT'
        end
        [UI.FD{n}.f,UI.FD{n}.val] = FFT_time2freq( t,val );
    else
        UI.FD{n}.f = freq;
        if strcmpi(AR_order,'auto')
            AR_order = 2;
            EC = -1;
            while 1
                [val_ar t_ar UI.FD{n}.val EC] = AR_estimate( t, val, freq, AR_order);
                if (EC==11)
                    AR_order = AR_order*2;
                else
                    break;
                end
            end
            if (EC~=0)
                warning('CSXCAD:ReadUI','AR estimation failed, skipping...')
                UI.FD{n}.val = DFT_time2freq( t, val, freq, SignalType );
            end
        elseif (AR_order<=0)
            UI.FD{n}.val = DFT_time2freq( t, val, freq, SignalType );
        else
            [val_ar t_ar UI.FD{n}.val EC] = AR_estimate( t, val, freq, AR_order);
            if (EC~=0)
                warning('CSXCAD:ReadUI','AR estimation failed, skipping...')
                UI.FD{n}.val = DFT_time2freq( t, val, freq, SignalType );
            end
        end
    end   
end
