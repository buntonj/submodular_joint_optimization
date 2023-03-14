function [x_primal,x_dual,dual_values,primal_values,gaps,added1,subopt,maxiterflag]  = minimize_submodular_projected_subgradient_descent_lovasz_polyak(F,param_F,opt_params)
% Compute  min_{x in [0,1]^p} f(x) using projected subgradient descent
% USING POLYAK'S RULE
%
% if sfm=1, perform submodular function minimization (and adapt gaps and
% function values accordingly)
%
%
% INPUT
% F,param_F: submodular function
% maxiter: maximum number of iterations
%
% OUTPUT
% x: argmin. If sfm=1, then outputs the subset
% values of cost function
% gaps: certified optimaly gaps at each iteration

if isfield(opt_params,'thresh')
    thresh = opt_params.thresh;
else
    thresh = 1e-6;
end

maxiter = opt_params.maxiter;

p = param_F.p;

% compute Lipschitz constant
FV = F(1:p);
Fsingletons = zeros(p,1);
FVsingletons = zeros(p,1);

for i=1:p
    Vmi = 1:p; Vmi(i)=[];
    FVsingletons(i) = FV-F(Vmi);
    Fsingletons(i) = F(i);
end
% B = sqrt( sum( max(FVsingletons.^2, Fsingletons.^2) ) );
B = sqrt( sum( (FVsingletons - Fsingletons).^2 ) );
D = sqrt(p);

w = rand(p,1);

% known in advance!
w(Fsingletons<=0)=0;
w(FVsingletons>=0)=1;

bestvalue = Inf;
s_ave = 0;
costs = [];

for iter =1:maxiter
    [s,Fvalues,order] = greedy_algo_submodular(w,F);
    s_ave = (s_ave * (iter-1) + s)/(iter);

    % compute function values
    [a,b] = min(Fvalues);
    if min(a,0) < bestvalue
        if a < 0
            % allow empty set to be optimal
            Aopt = order(1:b);
            bestvalue = a;
        else
            bestvalue = 0;
            Aopt = [];
        end
    end
    costs = [costs, bestvalue];
    added1(iter)  = s'*w;
    primal_values(iter) = a;
    dual_values(iter) = sum(min(s_ave,0));
    gaps(iter) = primal_values(iter) - dual_values(iter);
    if gaps(iter) <= thresh
        break
    end
    added2(iter) = added1(iter) - dual_values(iter);
    
    % projected gradient descent (Polyak's rule)
    w = w - ( s'*w - max(dual_values) ) / norm(s)^2 *  s;
    w = min(max(w,0),1);

end
maxiterflag = iter==maxiter;
subopt = gaps(end);
A = Aopt;
x_primal = Aopt;
x_dual = s_ave;


