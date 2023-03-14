%% CODE minimizing the given continuous submodular function F
% Uses Bach's provided Pair-wise FW algorithm
% param MUST hold:
% k - the discretization parameter for the region
% n - dimensionality of decision variable
% all other parameters necessary to evaluate F
% opt_params MUST hold:
% maxiter - maximum allowed iterations
% thresh - stopping duality gap criteria, default 1e-8
function [I_star,x_star,subopt,maxiterflag,gaps,costs] = bach_solve(F,param,opt_params)

maxiter = opt_params.maxiter;
k = param.k;
n = param.n;
if isfield(opt_params,'thresh')
    thresh = opt_params.thresh;
else
    thresh = 1e-6;
end
% random initialization
H0 = F(zeros(1,n),param);
w = greedy_algorithm(fliplr(cumsum(rand(n,k-1),2)),F,param);

ws = zeros( n*(k-1) ,maxiter+1 );
ws(:,1) = reshape(w,n*(k-1),1);
convex_combinations = zeros(1,maxiter+1);
convex_combinations(1) = 1;

%fprintf('Pairwise Frank Wolfe - 400 iterations\n');
costs = [];
for iter=1:maxiter
    % compute gradient direction
    rho = zeros(n,k-1);
    for i=1:n
        rho(i,:) = -pav(w(i,:));
    end
    [theta,~] = theta_minimizer(rho,F,param);
    y = (theta/(k-1)*2-1)';
    if strcmp(param.function_choice,'mr')
        costs = [costs, sparse_mr_cost(y,param)];
    elseif strcmp(param.function_choice,'int')
        costs = [costs, smooth_int_cost(y,param)];
    end
    
    % linear oracle
    [wbar,f,Fmin] = greedy_algorithm(rho,F,param);
    ws(:,iter+1) = reshape(wbar,n*(k-1),1);
    % compute away step
    ind = find(convex_combinations(1:iter)>0);
    [~,b] = min( ws(:,ind)'* reshape(rho,n*(k-1),1) );
    b = ind(b);
    away_direction = w - reshape(ws(:,b),n,k-1);
    max_step_away = convex_combinations(b);
    fw_direction = wbar - w;
    max_step_fw = 1;
    direction = fw_direction + away_direction;
    max_step = max_step_away;
    
    
    dual_fw_pair(iter) = .5 * sum( rho(:).^2 ) + sum( rho(:) .* w(:) );
    primal_fw_pair(iter) = f + .5 * sum( rho(:).^2 );
    primal_fw_pair_min(iter) = Fmin;
    dual_fw_pair_min(iter) = sum( min(min( cumsum(w,2) , [],2),0) );
    gaps(iter) = primal_fw_pair_min(iter)-dual_fw_pair_min(iter)-H0;
    if gaps(iter) < thresh
        break
    end
    % line search
    aa = sum( rho(:) .* ( direction(:) ) );
    bb = sum( ( direction(:) ).^2 );
    step = min(max_step,max(aa/bb,0));
    convex_combination_direction = zeros(1,maxiter+1);
    convex_combination_direction(iter+1)=1;
    convex_combination_direction(b) = -1;
    convex_combinations = convex_combinations + step * convex_combination_direction;
    w = w + step * direction;
    
end
maxiterflag = iter == maxiter;
subopt = primal_fw_pair_min(iter)-dual_fw_pair_min(iter)-H0;
[theta,~] = theta_minimizer(rho,F,param);
x_star = (theta/(k-1)*2-1)';
I_star = find(x_star);
end