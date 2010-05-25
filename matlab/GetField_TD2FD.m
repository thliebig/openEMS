function field_FD = GetField_TD2FD(field_TD, freq)
% function field_FD = GetField_TD2FD(field_TD, freq)
%
% transform time-domain field data into frequency domain
% 
% example:
%   freq = linspace(0,1e9,100); %target frequency vector
%   field = ReadHDF5FieldData('tmp/Ht.h5');
%   field_FD = GetField_TD2FD(field, freq);
%
% See also ReadHDF5FieldData
%
% openEMS matlab interface
% -----------------------
% author: Thorsten Liebig

t = field_TD.time;

field_FD.freq = freq;

for nf = 1:numel(freq)
    field_FD.values{nf} = 0;
end
    
numTS = numel(field_TD.values);

for n=1:numTS
    for nf = 1:numel(freq)
        f = freq(nf);
        field_FD.values{nf} = field_FD.values{nf} + 2/numTS * field_TD.values{n}.*exp(-1i*2*pi*f*t(n));
    end
end

