function F = modified_range(I,params)
% Modified range function + cardinality
if any(I)
    F = length(params.w) + max(I) - min(I) + length(sfo_unique_fast(I));
else
    F = 0;
end
end