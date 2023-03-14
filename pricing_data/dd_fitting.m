P = [avg_price; ones(1,num_intervals)]';
D = num_sales';
lambda = 0.01;
N = (num_products+1)*num_products; % problem dimension

cvx_begin quiet
    variable B(num_products+1,num_products)
    minimize(norm(P*B-D)+lambda*norm(B,'fro'))
    subject to
        %-B(1:num_products,:) == semidefinite(num_products);
        diag(B) <= 0;
        % add diagonal dominance constraint
        for i = 1:num_products
            -2*B(i,i) - sum(abs(B(i,:))) >= 0;
            -2*B(i,i) - sum(abs(B(1:num_products,i))) >= 0;
        end
cvx_end

beta = B(1:num_products,1:num_products); %pull out linear term of regression fit
alpha = B(end,1:num_products); % pull out bias of regression fit

max_prices = max(avg_price,[],2); % maximum price of each item
cost = 0.4; % percent of max price that each item costs
c = max_prices*0.4; % create costs vector
p_lower = 0*max_prices; % lower bound on retail price for each item, left zero here
% price i set to zero => no effect on other products' demand (since beta(i,j)*p(i) = 0) AND no revenue (NEGATIVE revenue actually!).

% building quadratic form for price optimization
% f(p) = p'*H*p + b'*p + const;
% want to minimize f with p >= p_lower.
Q = -0.5*(beta + beta'); % quadratic term in cost
b = -(2*p_lower'*Q + alpha + c'*beta)';
const = -c'*beta*p_lower + c*alpha;

%problem scaling
sc = norm(Q);
Q = Q/sc;
b = b/sc;