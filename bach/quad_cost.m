function val = quad_cost(y,param,set_fn)
val = 0.5*y*param.H*y' + y*param.b + param.lambda*set_fn(find(y));
end