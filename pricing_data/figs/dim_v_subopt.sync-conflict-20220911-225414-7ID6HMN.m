figure;
hold on;
for i = 3:3
    semilogy(dims,mnp_bounds(:29,i),'black','linewidth',2);
    semilogy(dims,(pgd_costs(:29,i)-mnp_costs(:,i)+mnp_bounds(:,i)),'red','linewidth',1);
    %semilogy(dims,mnp_costs(:,i)-mnp_bounds(:,i),'LineStyle',':','Color','#EDB120','linewidth',2)
end
iter_bound = ones(num_dims-1,1)*0 + eps;
semilogy(dims,iter_bound,'LineStyle',':','Color','#EDB120','linewidth',1.5)
legend({'Lift + MNP + FNNQP','PGD+FNNQP','Bound'});
title('Suboptimality over problem instances');
xlabel('Number of products (n)');
ylabel('Cost vs. Bound');
xlim([dims(1)-1,dims(end)+1]);
grid on;