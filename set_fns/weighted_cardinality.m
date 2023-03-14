function F = weighted_cardinality(I,params)
    F = sum(params.w(I));
end