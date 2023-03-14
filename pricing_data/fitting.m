% avg_price = num_products x num_intervals
% num_sales = num_products x num_sales

P = [avg_price; ones(1,num_intervals)]';
D = num_sales';
lambda = 1.0; % regularization for ridge regression.  higher reg -> "more submodular" quadratic form
B = pinv(P'*P + lambda*eye(num_products+1,num_products+1))*P'*D; %solve ridge regression problem

beta = B(1:num_products,1:num_products);
alpha = B(end,1:num_products);

max_prices = max(avg_price,[],2); % maximum price of each item
cost = 0.4; % percent of max price that each item costs
c = max_prices*0.4; % create costs vector
p_lower = cost*0.5*max_prices; % lower bound on price for each item (i.e. we will not take larger loss than half)

% building quadratic form for price optimization
H = -0.5*(beta + beta'); % quadratic term in cost
b = -(2*p_lower'*H + alpha + c'*beta);
const = -c'*beta*p_lower + c*alpha;

% COST FUNCTION IS:
% f(z) = z'*H*z + b'*z + const
% z >= 0