function FDTD = SetGaussExcite(FDTD,f0,fc,varargin)
% function FDTD = SetGaussExcite(FDTD,f0,fc,varargin);
%
% - f0 : center frequency
% - fc : 20dB cutoff frequency --> bandwidth is 2*fc
%
% optional key/value arguments:
% - 'DCFree', 1 : use a sine instead of a cosine carrier, so the pulse has
%                 no DC part and leaves no static charge behind. The spectrum
%                 falls to zero at DC below about f0-fc. With f0=0 the pulse
%                 is the derivative of the Gaussian, peak at about 0.47*fc.
%
% see also SetSinusExcite SetCustomExcite
%
% e.g
%
%     FDTD = SetGaussExcite(FDTD,1e9,1e8);
%     FDTD = SetGaussExcite(FDTD,1.55e9,1.45e9,'DCFree',1);
%
% openEMS matlab interface
% -----------------------
% author: Thorsten Liebig

FDTD.Excitation.ATTRIBUTE.Type=0;
FDTD.Excitation.ATTRIBUTE.f0=f0;
FDTD.Excitation.ATTRIBUTE.fc=fc;
FDTD.ATTRIBUTE.f_max=f0+fc;

for n=1:2:numel(varargin)
    if strcmpi(varargin{n},'DCFree')
        FDTD.Excitation.ATTRIBUTE.DCFree=double(varargin{n+1}~=0);
    else
        error('openEMS:SetGaussExcite',['unknown argument: ' varargin{n}]);
    end
end
