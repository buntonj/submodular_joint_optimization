function y = disc_to_cont_x(x,param)
n = length(x);
y = x / (param.k - 1) * 2 - 1;
end