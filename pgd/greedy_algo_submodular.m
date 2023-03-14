function [x,Fvalues,I] = greedy_algo_submodular(w,F);
% greedy algorithm started w
n = length(w);
x = zeros(n,1); [~,I] = sort(w,'descend');
Fvalues = zeros(n,1);
Fvalues(1) = F(I(1));
x(I(1)) = Fvalues(1);
for i=2:n
    Fvalues(i) = F(I(1:i));
    x(I(i)) = Fvalues(i) - Fvalues(i-1);
end
