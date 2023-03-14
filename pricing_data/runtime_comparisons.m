n_min = 30; % smallest problem size
n_max = 30; % largest problem size
step = 2; % step between instances
seed = RandStream('mlfg6331_64');

dims = n_min:step:n_max;
num_dims = length(dims);
num_runs = 4; % how many runs to do at each problem size?

pgd_times = zeros(num_dims,num_runs);
pgd_costs = zeros(num_dims,num_runs);
mnp_times = zeros(num_dims,num_runs);
mnp_costs = zeros(num_dims,num_runs);
mnp_bounds = zeros(num_dims,num_runs);

num_subs = zeros(num_dims,num_runs);
num_sups = zeros(num_dims,num_runs);

if exist('data')~= 1
    importing; % import the data if we haven't already
end

for i_d = 1:num_dims
    num_products = dims(i_d); % set the number of products to the current problem size
    fprintf('Running simulations for n = %d\n',num_products);
    for i_r = 1:num_runs
        fprintf('run %d',i_r);
        cleaning;
        fprintf('c');
        dd_fitting;
        fprintf('f\n');
        optimizing;
        % save timings
        pgd_times(i_d,i_r) = time_pgd;
        mnp_times(i_d,i_r) = time_mnp;
        % save costs
        pgd_costs(i_d,i_r) = cost_pgd;
        mnp_costs(i_d,i_r) = cost_mnp;
        %if cost_pgd > cost_mnp
        %    fprintf('!');
        %elseif cost_mnp == cost_pgd
        %    fprintf('=');
        %else
        %    fprintf('.');    
        %end
        %fprintf('\n');
        mnp_bounds(i_d,i_r) = bound;
        num_subs(i_d,i_r) = num_sub;
        num_sups(i_d,i_r) = num_sup;
    end
    %fprintf('\n');
end
pgd_subopt = pgd_costs-(mnp_costs-mnp_bounds);
mnp_subopt = mnp_bounds;

avg_pgd_times = mean(pgd_times,2);
avg_mnp_times = mean(mnp_times,2);
