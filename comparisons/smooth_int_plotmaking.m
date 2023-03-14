%% Plotting
true_signal = signal;
idxs = 1:length(signal);
noisy_signal = -param.b;
num_ns = length(dims);

%% Signal comparisons
% NOISY vs TRUE
sc = 1.1;
figure('Position',[0,0,750*sc,750*sc]);
subplot(4,1,1)
plot(idxs,true_signal,'linewidth',2);
hold on;
plot(idxs,noisy_signal,':','linewidth',2);
title('Provided Signal');
legend({'Noiseless Signal','Noisy Signal'});
xlabel('vector index');
ylabel('$\mathbf{y}$','Interpreter','latex');
xlim([1,n]);
grid on;

% MNP RECONSTRUCTION
subplot(4,1,2)
plot(idxs,true_signal,'linewidth',2);
hold on;
plot(idxs,x_star_mnp_fnnqp,'black','linewidth',2);
legend('Noiseless Signal','MNP + CPLEX/FNNQP');
title('Min Norm + CPLEX/FNNQP Reconstruction');
xlabel('vector index');
ylabel('$\mathbf{x}$','Interpreter','latex');
xlim([1,n]);
grid on;

% PGD RECONSTRUCTION
subplot(4,1,3)
plot(idxs,true_signal,'linewidth',2);
hold on;
plot(idxs,x_star_pgm_fnnqp,'red','linewidth',2);
title('Projected (Sub)Gradient Reconstruction');
legend('Noiseless Signal','PGD + CPLEX/FNNQP');
xlabel('vector index');
ylabel('$\mathbf{x}$','Interpreter','latex');
xlim([1,n]);
grid on;

% CONT SM RECONSTRUCTION
subplot(4,1,4)
plot(idxs,true_signal,'linewidth',2);
hold on;
plot(idxs,x_star_bach,'blue','linewidth',2);
title('Cont Submodular Reconstruction');
legend('Noiseless Signal','Cont Submodular');
ylabel('$\mathbf{x}$','Interpreter','latex');
xlabel('vector index');
xlim([1,n]);
grid on;
savefig('signal_comps_smooth_int.fig');


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
legend('PGD+CPLEX/FNNQP','Cont Submodular','MNP + CPLEX/FNNQP');
%set(gca,'fontsize',16);
%xlim([0,50]);
ylim([10^0.4,10^1.5]);
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
ylim([10^-6,10^5]);
xlim([8,112]);
xlabel('Problem Dimension (n)');
ylabel('Time (s)');
legend({'PGD + CPLEX','PGD+FNNQP','Cont Submodular','MNP + CPLEX','MNP+FNNQP'},'Location','southeast');
title('Runtime Scaling Comparison');

savefig('costs_comp_smooth_int.fig');