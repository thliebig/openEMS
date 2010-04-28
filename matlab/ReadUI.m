function UI = ReadUI(files, path)
% function UI = ReadUI(files, path)
%
% read current and voltages from multiple files found in path
%
% returns voltages/currents in time and frequency-domain
%
% e.g.
% UI = ReadUI({'ut1_1','ut1_2','it1'},'tmp/');
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
    tmp = load([path filenames{n}]);
    t = tmp(:,1)';
    val = tmp(:,2)';

    UI.TD{n}.t = t;
    UI.TD{n}.val = val;
    
    [UI.FD{n}.f,UI.FD{n}.val] = FFT_time2freq( t,val );
end
