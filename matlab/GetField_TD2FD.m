function field_FD = GetField_TD2FD(field_TD, freq)
% function field_FD = GetField_TD2FD(field_TD, freq)
%
% Transforms time-domain field data into the frequency domain
% Autocorrects the half-timestep offset of the H-field
% 
% example:
%   freq = linspace(0,1e9,100); %target frequency vector (Hz)
%   field = ReadHDF5FieldData('tmp/Ht.h5');
%   field_FD = GetField_TD2FD(field, freq);
%
% openEMS matlab interface
% -----------------------
% author: Thorsten Liebig
%
% See also ReadHDF5FieldData

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
        % t(n) is absolute time and therefore the half-timestep offset of
        % the H-field is automatically compensated
        % openEMS output: E-fields start at t=0
        % openEMS output: H-fields start at t=delta_t/2
    end
end

