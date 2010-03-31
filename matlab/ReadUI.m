function UI = ReadUI(files, path)

if (nargin<2)
    path ='';
end

if (ischar(files))
    filenames{1}=files;
else
    filenames=files;
end

for n=1:numel(filenames)
    tmp = load([path filenames{n}]);
    t = tmp(:,1)';
    val = tmp(:,2)';

    UI.TD{n}.t = t;
    UI.TD{n}.val = val;
    
    [UI.FD{n}.f,UI.FD{n}.val] = FFT_time2freq( t,val );
end
