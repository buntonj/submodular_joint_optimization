num_products = 20;
%importing;
cleaning;
dd_fitting;
optimizing;

num_mnp_iters = length(costs_mnp);
new_costs_mnp = zeros(num_mnp_iters,1);
for i = 1:num_mnp_iters
    Icurr = lifted_prob.order_reversing_bijection(find(Is(i,:)));
    Iz = Icurr(Icurr <= n);
    Iw = Icurr(Icurr > n) - n;
    Iu = intersect(Iz,Iw);
    Iv = union(Iz,Iw);
    cz = prob.joint_value(Iz);
    cw = prob.joint_value(Iw);
    cu = prob.joint_value(Iu);
    cv = prob.joint_value(Iv);
    new_costs_mnp(i) = min([cz,cw,cu,cv]);
end

cost_lower_bound = cost_mnp - bound;
subopt_mnp = new_costs_mnp- cost_lower_bound + eps;
subopt_mnp = [subopt_mnp; subopt_mnp(end)*ones(100-length(subopt_mnp),1)];
subopt_pgd = costs_pgd-cost_lower_bound + 100*eps;
iter_bound = ones(100,1)*bound + eps;
figure;
semilogy(subopt_mnp(1:100),'black','linewidth',2);
hold on;
semilogy(subopt_pgd,'red','linewidth',2);
semilogy(iter_bound,':','linewidth',2);
legend({'Lift + MNP + FNNQP','PGD + FNNQP','Bound'});
title('Suboptimality over Iterations');
xlabel('Iteration');
ylabel('Suboptimality');
grid on;