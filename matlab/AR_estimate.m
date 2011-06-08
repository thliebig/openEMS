function [val_ar t_ar f_val_ar] = AR_estimate( t, val, freq, nu, mu, expand_factor)
% [val_ar t_ar f_val_ar] = AR_estimate( t, val, freq, < nu, mu, expand_factor >)
%
% apply autoregressive signal model to improve dft results
%  
% t     : time vector
% val   : time domain values
% freq  : frequency vector for dft
%
% optional
% nu    : AR order (default 40)
% mu    : number of timesteps to train the model (default 3*nu)
% expand_factor : increase signal length by this factor (default 5)
%
% 
% openEMS matlab interface
% -----------------------
% Author: Thorsten Liebig, 2011
%
% See also ReadUI, DFT_time2freq

if numel(t) ~= numel(val)
    error 'numel(t) ~= numel(val)'
end

if (nargin<4)
    nu = 40;
end
if (nargin<5)
    mu = 3*nu;
end
if (nargin<6)
    expand_factor=5;
end

if (mu<=2*nu)
    error 'mu should be larger than 2*nu'
end

if (expand_factor<=1)
    error 'expand_factor must be larger than 1'
end

dt = t(2)-t(1);

M = numel(t);

if (M<0.6*mu)
    error 'signal is to short for AR estimate --> decrease AR order'
end

for n=1:mu-nu
    b(n) = val(end-n+1);
    for m=1:nu
        A(n,m)=val(end-n+1-m);
    end
end

a = ((A'*A)\A')*b';

val_ar = val;
t_ar = t;
for k=M:expand_factor*M
    val_ar(k) = 0;
    t_ar(k) = t_ar(k-1)+dt;
    val_ar(k) = sum(a.*val_ar(k-(1:nu))');
end

if (max(val_ar(M:end)) > max(val))
    error 'estimated signal appears to be unstable --> use a larger mu'
end

f_val_ar = DFT_time2freq(t_ar, val_ar, freq);
