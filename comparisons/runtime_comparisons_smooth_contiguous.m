%% initialization
clear

addpath('../problems'); % joint problem object directory
addpath('../pgd'); % projected subgradient descent algorithm directory
addpath('../bach'); % continuous sm min algorithm directory
addpath(genpath('../set_fns')); % set function regularizers directory

seed = 2;
rng(2,'twister');

% dimensions to consider
dim_start = 10; %starting dimension
dim_step = 2; %step size between dimensions
num_steps = 50; %number of steps up in dimension
dims = dim_start:dim_step:(dim_start + num_steps*dim_step);

%% Stats to track
%times = zeros(num_steps+1,5); % runtimes for each method
bach_times = zeros(num_steps+1,1);
mnp_cplex_times = zeros(num_steps+1,1);
mnp_fnnqp_times = zeros(num_steps+1,1);
pgm_cplex_times = zeros(num_steps+1,1);
pgm_fnnqp_times = zeros(num_steps+1,1);

costs = zeros(num_steps+1,5); % final incurred costs by each solution
disagree = zeros(num_steps+1,2); % indices where continuous or pgm algorithm disagree with us

%final suboptimality gaps
mnp_cplex_subopt = zeros(num_steps+1,1);
mnp_fnnqp_subopt = zeros(num_steps+1,1);
bach_subopt = zeros(num_steps+1,1);
pgm_cplex_subopt = zeros(num_steps+1,1);
pgm_fnnqp_subopt = zeros(num_steps+1,1);

%flagging if continuous/pgm reached max iterations before converging
bach_maxiterflag = false(num_steps+1,1);
pgm_cplex_maxiterflag = false(num_steps+1,1);
pgm_fnnqp_maxiterflax = false(num_steps+1,1);
param.function_choice = 'int';

%% Iterate over dimensions
for dim = 1:num_steps+1
    %% Set up problem statement
    n = dims(dim); % dimension
    param.n = n;

    fprintf('Setting up problem parameters, N = %g\n',n);
    D = spdiags([zeros(n-1,1),ones(n-1,1),-ones(n-1,1)],-1:1,n-1,n); %differencing matrix
    mu = 0.8;
    param.H = 2*mu*(D'*D) + eye(n);
    t = (0:(n-1))/(n-1)*2-1;
    signal = 100 *(  (t<-.2) .* (t>-.8) .* ((t+.5).^2-.3^2).^2 + (t>.2) .* (t<.8) .* ((t-.5).^2-.3^2).^2 );
    param.b =  -(signal + 0.1*randn(1,n))';
    param.lambda = 0.05;
    param.wt = ones(1,n);

    opt_params.maxiter = 100; %maximum number of iterations
    opt_params.thresh = 1e-6;


    %% Frank-wolfe on B(F) with pairwise steps
    % from Lacoste-Julien and Jaggi (2015)
    fn_params = struct();
    fn_params.w = param.wt;
    fn_params.choice = param.function_choice;
    set_fn = num_intervals_handle(fn_params);
    %continuous_F = @continuous_submodular_function_mr;
    continuous_F = @continuous_submodular_function_int;
    %continuous_F = @(x,param)( quad_cost(disc_to_cont_x(x,param),param,set_fn));
    param.k = 51; % elements/variable (discretization of distributions on each dimension)
    fprintf('Running CONT-SFM algorithm...\n');
    tic
    [I_star_bach,x_star_bach,bach_subopt(dim),bach_maxiterflag(dim),subopt_path_bach,costs_bach] = bach_solve(continuous_F,param,opt_params);
    bach_times(dim) = toc;

    
    %% Our Method using the MNP algorithm and CPLEX
    fprintf('Running MNP + CPLEX algorithm...\n');
    fn_params = struct();
    fn_params.w = param.wt;
    fn_params.choice = param.function_choice;
    set_fn = num_intervals_handle(fn_params);
    prob = cplex_joint_problem(param.H,param.b,param.lambda,set_fn);
    prob.opt = sfo_opt({'verbosity_level',0,'minnorm_stopping_thresh',opt_params.thresh});
    tic
    prob.prune_space();
    [I_star_mnp_cplex,x_star_mnp_cplex,subopt_path_mnp_cplex,costs_mnp_cplex] = prob.optimize();
    mnp_cplex_subopt(dim) = subopt_path_mnp_cplex(end);
    mnp_cplex_times(dim) = toc;

    %% Our method using MNP algorithm and FNNQP for sub-problems
    fprintf('Running MNP + FNNQP algorithm...\n');
    prob2 = fnnqp_joint_problem(param.H,param.b,param.lambda,set_fn);
    prob2.opt = sfo_opt({'verbosity_level',0,'minnorm_stopping_thresh',opt_params.thresh});
    tic
    prob2.prune_space();
    [I_star_mnp_fnnqp,x_star_mnp_fnnqp,subopt_path_mnp_fnnqp,costs_mnp_fnnqp] = prob2.optimize();
    mnp_fnnqp_subopt(dim) = subopt_path_mnp_fnnqp(end);
    mnp_fnnqp_times(dim) = toc;
    
    %% Projected Subgradient Descent with Polyak's Rule, CPLEX
    % Runs on the parameterized function, still sort of "Our Method"
    
    fprintf('Running PGD + CPLEX...\n');
    F_joint = sfo_fn_wrapper(@(A) prob.joint_value(A));
    F_param.p = n;
    prob.bottom = [];
    prob.V = 1:n;
    tic
    prob.prune_space();
    [I_star_pgm_cplex,~,~,costs_pgm_cplex,subopt_path_pgm_cplex,~,pgm_subopt_cplex(dim),pgm_maxiterflag_cplex(dim)] = minimize_submodular_projected_subgradient_descent_lovasz_polyak(F_joint,F_param,opt_params);
    x_star_pgm_cplex = prob.set_to_minimizer(I_star_pgm_cplex);
    %}
    pgm_cplex_times(dim) = toc;
    
    %% Projected Subgradient Descent with Polyak's Rule, FNNQP
    fprintf('Running PGD + FNNQP...\n');
    F_joint = sfo_fn_wrapper(@(A) prob2.joint_value(A));
    F_param.p = n;
    prob2.bottom = [];
    prob2.V = 1:n;
    tic
    prob2.prune_space();
    [I_star_pgm_fnnqp,~,~,costs_pgm_fnnqp,subopt_path_pgm_fnnqp,~,pgm_subopt_fnnqp(dim),pgm_maxiterflag_fnnqp(dim)] = minimize_submodular_projected_subgradient_descent_lovasz_polyak(F_joint,F_param,opt_params);
    x_star_pgm_fnnqp = prob2.set_to_minimizer(I_star_pgm_fnnqp);
    %}
    pgm_fnnqp_times(dim) = toc;
    
    disagree(dim,1) = sum((x_star_bach~=0) ~= (x_star_mnp_cplex~=0));
    disagree(dim,2) = sum((x_star_pgm_cplex~=0) ~= (x_star_mnp_cplex~=0));
    costs(dim,:) = [prob2.joint_value(I_star_mnp_cplex), prob2.joint_value(I_star_mnp_fnnqp), prob2.joint_value(I_star_bach), prob2.joint_value(I_star_pgm_cplex), prob2.joint_value(I_star_pgm_fnnqp)];
    %times(dim,:) =[mnp_cplex_time, mnp_fnnqp_time, bach_time, pgm_cplex_time, pgm_fnnqp_time];
end
costs_mnp_cplex = costs_mnp_cplex + 0.5*param.b'*param.b;
costs_mnp_fnnqp = costs_mnp_fnnqp + 0.5*param.b'*param.b;
costs_bach = costs_bach + 0.5*param.b'*param.b;
costs_pgm_cplex = costs_pgm_cplex + 0.5*param.b'*param.b;
costs_pgm_fnnqp = costs_pgm_fnnqp + 0.5*param.b'*param.b;