function field = GetField_TD2FD(field, freq)
% function field = GetField_TD2FD(field, freq)
%
% Transforms time-domain field data into the frequency domain
% Auto-corrects the half-timestep offset of the H-field
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

if (~isfield(field,'TD'))
    warning('openEMS:GetField_TD2FD','field has no time domain data... skipping FD transformation...');
    return
end

t = field.TD.time;
dt = t(2)-t(1);

clear field.FD

field.FD.frequency = freq;

for nf = 1:numel(freq)
    field.FD.values{nf} = 0;
end
    
numTS = numel(field.TD.values);

for n=1:numTS
    for nf = 1:numel(freq)
        f = freq(nf);
        field.FD.values{nf} = field.FD.values{nf} + field.TD.values{n}.*exp(-1i*2*pi*f*t(n)) * 2 * dt;
        % t(n) is absolute time and therefore the half-timestep offset of
        % the H-field is automatically compensated
        % openEMS output: E-fields start at t=0
        % openEMS output: H-fields start at t=delta_t/2
    end
end

field.FD.DataType=1;

