function f = continuous_submodular_function_mr(x,param)
% 1D Laplacian and concave function to induce sparsity
% parameters
% k: cardinality of each variable
% lambda: regularization parameter for sparsity
% mu : regularization parameter for smoothness
n = length(x);
y = x / (param.k - 1) * 2 - 1;
%y = x'; % for debugging
% f = sum( ( y - param.signal ).^2 ) + param.mu * sum( (y(2:n)-y(1:n-1)).^2) + param.lambda * sum(min(abs(y),1/(param.k-1)));
%f = sum( ( y - param.signal ).^2 ) + param.mu * sum( (y(2:n)-y(1:n-1)).^2) + param.lambda * sum(abs(y).^.125);

supp = find(y);
% Modified Range Function + cardinality
if any(supp)
    submod_part = (length(supp) + max(supp) - min(supp) + length(y));
else
    submod_part = 0;
end

% Weighted cardinality
% submod_part = sum(param.wt(y~=0));

% Number of intervals + cardinality
%d = (supp-circshift(supp,-1)) ~= -1;
%submod_part = sum(param.wt(supp)) + sum(d(1:(length(d)-1))) + 1;

f = 0.5*y*param.H*y' + y*param.b + param.lambda*submod_part;
end

