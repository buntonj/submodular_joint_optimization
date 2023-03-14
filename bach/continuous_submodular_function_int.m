function f = continuous_submodular_function_int(x,param)
% parameters
% k: cardinality of each variable
% lambda: regularization parameter for sparsity
% mu : regularization parameter for smoothness
n = length(x);
y = x / (param.k - 1) * 2 - 1;
%y = x';
% f = sum( ( y - param.signal ).^2 ) + param.mu * sum( (y(2:n)-y(1:n-1)).^2) + param.lambda * sum(min(abs(y),1/(param.k-1)));
%f = sum( ( y - param.signal ).^2 ) + param.mu * sum( (y(2:n)-y(1:n-1)).^2) + param.lambda * sum(abs(y).^.125);

supp = find(y);
% Modified Range Function + cardinality
%if isempty(supp)
%    submod_part = 0;
%else
%    submod_part = sum(param.wt(y~=0)) + max(supp) - min(supp) + length(y);
%end

% Weighted cardinality
% submod_part = sum(param.wt(y~=0));

% Number of intervals + cardinality
d = (supp-circshift(supp,-1)) ~= -1;
submod_part = sum(param.wt(supp)) + sum(d(1:(length(d)-1))) + 1;

%d = (supp - circshift(supp,-1)) ~= 0; %1 when there are successive entries zero vs nonzero
%f = 0.5*y*param.H*y' + y*param.b + param.lambda*(sum(param.wt(supp)) + sum(d(1:(length(d)-1))));
%f = 0.5*y*param.H*y' + y*param.b + param.lambda*(sum(param.wt(supp)));
f = 0.5*y*param.H*y' + y*param.b + param.lambda*submod_part;
end

