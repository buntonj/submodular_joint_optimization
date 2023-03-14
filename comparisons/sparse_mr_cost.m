function cost = sparse_mr_cost(x,param)
continuous_part  = 0.5*x'*param.H*x + param.b'*x;
supp = find(x);
submodular_part = any(supp)*(max(supp)-min(supp) + length(x) + length(supp));
cost = continuous_part + param.lambda*submodular_part;
end