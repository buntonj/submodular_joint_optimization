% making plots

true_signal = signal;
mnp_cplex_signal = param.R*x_star_mnp_cplex;
mnp_fnnqp_signal = param.R*x_star_mnp_fnnqp;
pgd_cplex_signal = param.R*x_star_pgm_cplex;
pgd_fnnqp_signal = param.R*x_star_pgm_fnnqp;
bach_signal = param.R*x_star_bach;

sc = 1.1;
figure('Position',[0,0,750*sc,750*sc]);

subplot(4,1,1)
plot(true_signal,'linewidth',2);
title('Provided Signal');
legend('Signal');
xlabel('vector index');
ylabel('$\mathbf{b}$','Interpreter','latex');
ylim([0,1.1]);
xlim([1,n]);
grid on;

subplot(4,1,2)
plot(true_signal,'linewidth',2);
hold on;
plot(mnp_fnnqp_signal,'black','linewidth',2);
legend('Signal','MNP + CPLEX/FNNQP');
title('Min Norm + CPLEX Reconstruction');
xlabel('vector index');
ylabel('$\mathbf{Dx}$','Interpreter','latex');
ylim([0,1.1]);
xlim([1,n]);
grid on;

subplot(4,1,3)
plot(true_signal,'linewidth',2);
hold on;
plot(pgd_fnnqp_signal,'red','linewidth',2);
legend('Signal','PGD + CPLEX/FNNQP');
title('Projected (Sub)Gradient Reconstruction');
xlabel('vector index');
ylabel('$\mathbf{Dx}$','Interpreter','latex');
xlim([1,n]);
ylim([0,1.1]);
grid on;

subplot(4,1,4)
plot(true_signal,'linewidth',2);
hold on;
plot(bach_signal,'blue','linewidth',2);
legend('Signal','Cont Submodular');
title('Cont Submodular Reconstruction');
xlabel('vector index');
ylabel('$\mathbf{Dx}$','Interpreter','latex');
xlim([1,n]);
ylim([0,1.1]);
grid on;

%% COSTS PLOT
n_costs = opt_params.maxiter;
costs_mnp_cplex_plots = [costs_mnp_cplex, costs_mnp_cplex(end)*ones(1,n_costs-length(costs_mnp_cplex))];
costs_mnp_fnnqp_plots = [costs_mnp_fnnqp, costs_mnp_fnnqp(end)*ones(1,n_costs-length(costs_mnp_fnnqp))];
costs_bach_plots = [costs_bach, costs_bach(end)*ones(1,n_costs-length(costs_bach))];
costs_pgm_cplex_plots = [costs_pgm_cplex, costs_pgm_cplex(end)*ones(1,n_costs-length(costs_pgm_cplex))];
costs_pgm_fnnqp_plots = [costs_pgm_fnnqp, costs_pgm_fnnqp(end)*ones(1,n_costs-length(costs_pgm_fnnqp))];
figure('Position',[0,0,750*sc,250*sc]);
subplot(1,2,1)
semilogy(1:n_costs,max(costs_pgm_fnnqp_plots,eps),'red','linewidth',2);
hold on;
semilogy(1:n_costs,max(costs_bach_plots,eps),'blue','linewidth',2);
semilogy(1:n_costs,max(costs_mnp_fnnqp_plots,eps),'black','linewidth',1.5);
legend('PGD + CPLEX/FNNQP','Cont Submodular','MNP + CPLEX/FNNQP');
%set(gca,'fontsize',16);
%xlim([0,50]);
ylim([10^3.0301,10^3.06]);
xlabel('Iteration');
ylabel('Cost');
title('Objective Value');
grid on;

%% RUNTIMES PLOT
subplot(1,2,2)
semilogy(dims,pgm_cplex_times,'red','linewidth',2);
hold on;
grid on;
semilogy(dims,pgm_fnnqp_times,'LineStyle',':','Color','red','linewidth',2);
semilogy(dims,bach_times,'blue','linewidth',2)
semilogy(dims,mnp_cplex_times,'black','linewidth',2)
semilogy(dims,mnp_fnnqp_times,'LineStyle',':','Color','black','linewidth',2);
ylim([10^-6.5,10^5]);
xlim([8,112]);
xlabel('Problem Dimension (n)');
ylabel('Time (s)');
legend({'PGD + CPLEX','PGD + FNNQP','Cont Submodular','MNP + CPLEX','MNP+FNNQP'},'Location','southeast');
title('Runtime Scaling Comparison');