function [w,f,Fmin] = greedy_algorithm(rho,F,param_F)
% greedy algorithm (assumes all variables have the same cardinality)

[ n k ] = size(rho);
w = zeros(size(rho));
k = k + 1;

% first order all rhos (does preserve the ordering within rows if equal
% values)
[~, s] = sort( reshape(rho',n*(k-1),1), 'descend');
[js,is] = ind2sub([k-1 n],s);

% now go through all elements
xold = zeros(1,n);
Fold = F(xold,param_F);
H0 = Fold;

xmin = xold;
Fmin = Fold;


for i=1:n*(k-1)
    xnew = xold; xnew(is(i))=js(i);
    Fnew = F(xnew,param_F);
    if (Fnew<Fmin), xmin = xnew; Fmin = Fnew; end
    w(is(i),js(i)) = Fnew - Fold;
    xold = xnew;
    Fold = Fnew;
end
f = sum( w(:) .* rho(:) )+H0;
