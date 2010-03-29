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
    
    dt=t(2)-t(1);
    val = [val zeros(1,5000)];
    L=numel(val);
    UI.FD{n}.f = (0:L-1)/L/dt;
    UI.FD{n}.f = UI.FD{n}.f(1:floor(L/2));
    UI.FD{n}.val = fft(val)/L;
    UI.FD{n}.val = UI.FD{n}.val(1:floor(L/2));
end