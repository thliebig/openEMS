function pass = check_frequency( f, val, f_upper, f_lower, rel_amplitude, type )

pass = true;
max1 = max(val);

if numel(f_upper) ~= numel(f_lower)
    error 'inconsistant vectors'
end

for n=1:numel(f_upper)
    f1 = f_lower(n);
    f2 = f_upper(n);
    f1_idx = interp1( f, 1:numel(f), f1, 'nearest' );
%    if f(f1_idx) < f1, f1_idx = f1_idx + 1; end
    f2_idx = interp1( f, 1:numel(f), f2, 'nearest' );
%    if f(f2_idx) > f2, f2_idx = f2_idx - 1; end

    if strcmp( type, 'inside' )
        if max( val(f1_idx:f2_idx) ) < max1 * rel_amplitude
            pass = false;
            return
        end
    elseif strcmp( type, 'outside' )
        if max( val(f1_idx:f2_idx) ) > max1 * rel_amplitude
            pass = false;
            return
        end
    else
        error 'unsupported operation'
    end
end
