%% Using data thus far, build optimization problem and solve it
% Incoming releveant info is
% f(p) = p'*Q'p + b'*p + const
% n = num_products = size(Q,1)
% p >= p_lower (set to zero by default)

addpath('../problems'); % joint_problems objects and MNP
addpath('../pgd'); % projected subgradient descent algorithms
addpath(genpath('../set_fns')); % set function regularizers

n = num_products;
print_optimality = true; % whether or not to print if we do better
return_sets = true;

%% building regularizer
num_groups = floor(n/1) ;%floor(num_products/min(5,num_products)); % how many groups will there be?
% to track groups, build num_groups x n mask
G = false(n, num_groups);

% pick a random group for each product
for i = 1:n
    G(i,randi(num_groups)) = true;
end

% determine group start-up costs
w = zeros(num_groups,1);
bias = 0.75;
for i = 1:num_groups
    w(i) = 0.35*sum(c(G(:,i))) + bias; % group start up some frac of combined costs
end

% to compute the cost of a set of indices S, we would do:
% g(S) =  sum(w(any(G(S,:),1)));
% any(G(S,:),1) finds group indices of groups whose indices are in S
% summing w over these gives the total startup costs
fn_params = struct();
fn_params.w = w;
fn_params.G = G;
fn_params.choice = 'groups';
lambda = 1.0*sqrt(n/20);
set_fn = groups_handle(fn_params);

%% Constructing original joint problem
prob = fnnqp_joint_problem(Q,b,lambda,set_fn); % builds a joint problem instance

%% Constructing lifted joint problem
% Performing the lift!

% building lifted Q matrix
Qsub = diag(diag(Q)) + Q.*(Q <= 0);
Qsup = Q - Qsub;
lifted_Q = [Qsub, Qsup; Qsup, Qsub];
num_sup = sum(Qsup > 0,'all');
num_sub = sum(Qsub ~= 0,'all');

% lifted b matrix and associated constant terms
lifted_b = [b; b];
lifted_const = 2*const;

% checking that lifted problem is convex
check = min(eig(lifted_Q)); % if check < 0 then we can't use the lift
if (check < 0)
    fprintf('ERROR: Lifted quadratic is not convex.\n');
end
% create lifted problem object
lifted_prob = lifted_joint_problem(lifted_Q,lifted_b,lambda,set_fn);

%% Solve lifted problem with the MNP algorithm
%lifted_prob.prune_space();
tic
if return_sets
    [I_star_mnp,x_star_mnp,Is,costs_mnp] = lifted_prob.optimize_return_sets();
else
    [I_star_mnp,x_star_mnp,subopt,costs_mnp] = lifted_prob.optimize();
end
time_mnp =  toc;
flipped_I_star = lifted_prob.order_reversing_bijection(I_star_mnp);
I_z = flipped_I_star(flipped_I_star <= n);
I_w = flipped_I_star(flipped_I_star > n) - n;
z = x_star_mnp(1:n);
w = x_star_mnp(n+1:end);
z2 = min(z,w);
w2 = max(z,w);

% extract the best possible bound
bound = 0.5*(z-w)'*Qsup*(z-w);
bound = min(0.5*(z2-w2)'*Qsup*(z2-w2),bound);
bound = max(bound,0); % if by miracle or numerical error we are lower than zero, we're optimal
cost_mnp = min([prob.joint_value(I_z),prob.joint_value(I_w),prob.joint_value(find(z2)),prob.joint_value(find(w2))]);

%% Solve original joint problem with PGD
F_joint = sfo_fn_wrapper(@(A) prob.joint_value(A));
L = -prob.joint_value(prob.V) + 2*prob.submod_function(prob.V);
opt_params = struct();
opt_params.L = L;
opt_params.maxiter = 100;
V = 1:n;
tic;
[x_star_pgd,I_star_pgd,costs_pgd] = vanilla_pgd(F_joint,V,opt_params);
time_pgd = toc;
cost_pgd = costs_pgd(end);

if print_optimality
    if cost_mnp < cost_pgd ; 
        fprintf('MNP < PGD\n');
    elseif cost_mnp > cost_pgd
        fprintf('PGD < MNP\n');
    else
        fprintf('PGD = MNP\n');
    end
    fprintf('bound : %f \n',bound)
    fprintf('%d entries supermodular, ',num_sup)
    fprintf('largest supermodular entry %f  \n',max(Qsup,[],'all'))
    
end