function cost = smooth_int_cost(x,param)
continuous_part  = 0.5*x'*param.H*x + param.b'*x;
supp = find(x);
d = (supp-circshift(supp,-1)) ~= -1;
submodular_part = sum(param.wt(supp)) + sum(d(1:(length(d)-1))) + 1;
cost = continuous_part + param.lambda*submodular_part;
end