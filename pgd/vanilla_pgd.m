function [s, I_star, costs] = vanilla_pgd(F,V,opt_params)
T = opt_params.maxiter;
n  = length(V);
R = 2*sqrt(n);
L = opt_params.L;
eta = R/(L*sqrt(T));

costs = zeros(T,1);

s = zeros(1,n);

for i = 1:T
    w = sfo_polyhedrongreedy(F,V,s);
    costs(i) = sum(w.*s);
    s = s - eta*w; % gradient step
    s = min(max(0,s),1); % projection
end

[~,I] = sort(s,'descend');
best = F([]);
I_star = [];
for i = 1:n
    if F(I(1:i)) < best
        I_star = I(1:i);
        best = F(I(1:i));
    end
end
I_star = sort(I_star,'ascend');
costs(T) = best;

end