%used to check discretization dependence on a problem instance
k_dep = []; k_time = [];
i = 1;
for k = 51:20:401
param.k = k;
fprintf('Running for k = %g\n',k);
tic
[I_star_bach,x_star_bach,bach_subopt(dim),bach_maxiterflag(dim),subopt_path_bach,costs_bach] = bach_solve(continuous_F,param,opt_params);
k_dep(i) = costs_bach(end);
k_time(i) = toc;
i = i + 1;
end
k_dep = k_dep + 0.5*param.b'*param.b;
figure();
semilogy(51:20:401,k_dep,'blue','linewidth',3);
hold on;
semilogy(51:20:401,costs_us(end)*ones(size(k_dep)),'black','linewidth',3);
grid on;
xlabel('Discretization Level (k)');
ylabel('Cost');
title('Objective Value');
legend('Cont Submodular','Optimal')
set(gca,'fontsize',16)

figure();
semilogy(51:20:401,k_time,'blue','linewidth',3);
xlabel('Discretization Level (k)');
ylabel('Running Time (s)');
title('Runtime Discretization Dependence');
legend('Cont Submodular');
grid on;
set(gca,'fontsize',16)