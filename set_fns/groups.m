% GROUPS FUNCTION
function F = groups(I,params)
    F = sum(params.w(any(params.G(I,:),1)));
end