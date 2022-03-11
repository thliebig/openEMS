function [val_ar t_ar f_val_ar EC] = AR_estimate( t, val, freq, nu, mu, expand_factor)
% [val_ar t_ar f_val_ar  EC] = AR_estimate( t, val, freq, < nu, mu, expand_factor >)
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
% return values:
% val_ar:   AR estimated time signal
% t_ar:     time vector
% f_val_ar: FD transformed AR estimated signal
% EC:       error code
%           0 --> no error
%           1 --> input error: t and val mismatch
%           2 --> input error: mu has to be larger than 2*nu
%           3 --> input error: expand_factor has to be larger than 1
%           10 --> AR error: signal is to short for AR estimate --> decrease AR order
%           11 --> AR error: estimated signal appears to be unstable --> use a different mu
%
% openEMS matlab interface
% -----------------------
% Author: Thorsten Liebig, 2011
%
% See also ReadUI, DFT_time2freq

EC = 0;
val_ar = [];
t_ar = [];
f_val_ar = [];


if numel(t) ~= numel(val)
    if (nargout<4)
        error 'numel(t) ~= numel(val)'
    else
        EC = 1;
        return
    end
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
    if (nargout<4)
        error 'mu has to be larger than 2*nu'
    else
        EC = 2;
        return
    end
end

if (expand_factor<=1)
    if (nargout<4)
        error 'expand_factor has to be larger than 1'
    else
        EC = 3;
        return
    end
end

dt = t(2)-t(1);

M = numel(t);

if (M<0.6*mu)
    if (nargout<4)
        error 'signal is to short for AR estimate --> decrease AR order'
    else
        EC = 10;
        return
    end
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
    if (nargout<4)
        error 'estimated signal appears to be unstable --> use a different mu'
    else
        EC = 11;
        return
    end
end

f_val_ar = DFT_time2freq(t_ar, val_ar, freq);
