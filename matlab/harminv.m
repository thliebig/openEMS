function [f,decay,Q,amp,phase,err]=harminv( timeseries, timestep, start_freq, stop_freq )
% [f,decay,Q,amp,phase,err]=harminv( timeseries, timestep, start_freq, stop_freq )
%
% reconstruct time signal with:
% real(amp(n) * exp(-1i*(2*pi*freq(n)*t-phase(n))-decay(n)*t)));
%
% example:
% t = linspace(0,(1/0.3e9)*5,1000);
% u =     0.8 * sin(2*pi*0.7e9 * t + 0/180*pi);
% u = u + 0.3 * sin(2*pi*0.3e9 * t + 0/180*pi);
% [freq,decay,Q,amp,phase,err]=harminv( u, t(2)-t(1), 0, 1e9 );
%
% -----------------------
% Sebastian Held <sebastian.held@gmx.de>


%harminv_exe = '../../../harminv/harminv';
harminv_exe = 'harminv';
options = ['-t ' num2str(timestep)];
tmpfile = tempname;

% convert to column vector
if size(timeseries,2) > size(timeseries,1)
    timeseries = timeseries.';
end

% harminv hangs if timeseries is zero only
if all(timeseries == 0)
    disp( 'timeseries is 0' );
    return
end

% write timeseries to temporary file
dlmwrite( tmpfile, timeseries );

command = [harminv_exe ' ' options ' '...
           num2str(start_freq) '-' num2str(stop_freq) ' < ' tmpfile];
[status,result] = system( command );
if status ~= 0
    disp( 'error executing harminv:' );
    disp( command );
    disp( ['exit status: ' num2str(status)] );
    disp( ['output: ' result] );
    return
end

%disp( command )
%disp( result )

f = [];
decay = [];
Q = [];
amp = [];
phase = [];
err = [];

pos = 1;
while pos < length(result)
    result = result(pos:end);
    [C,pos] = textscan( result, '%f', 6, 'Delimiter', ',', 'HeaderLines', 1 );
    if isempty(C{1}), break; end
    if C{1}(1) >= 0
        f     = [f C{1}(1)];
        decay = [decay C{1}(2)];
        Q     = [Q C{1}(3)];
        amp   = [amp 2*C{1}(4)]; % neglecting negative frequencies => amplitude doubles
        phase = [phase C{1}(5)];
        err   = [err C{1}(6)];
    end
end
