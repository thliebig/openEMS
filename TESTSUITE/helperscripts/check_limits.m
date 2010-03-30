function pass = check_limits( Z, upper_limit, lower_limit )

% make row vector
if size(Z,1) ~= 1
    Z = Z.';
end

if numel(upper_limit) == 1
    upper_limit = upper_limit * ones(1,size(Z,2));
end
if numel(lower_limit) == 1
    lower_limit = lower_limit * ones(1,size(Z,2));
end


pass = 1;
if any( Z > upper_limit )
    pass = 0;
end
if any( Z < lower_limit )
    pass = 0;
end
