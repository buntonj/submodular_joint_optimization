percent_subs = 100*num_subs./(num_subs+num_sups);
cost_diffs = pgd_costs-mnp_costs;
figure;
hold on;
for i = 1:num_runs
    scatter(percent_subs(:,i),cost_diffs(:,i),'red')
end

plot([min(percent_subs,[],'all'),max(percent_subs,[],'all')],[0,0],'black','linewidth',2)
grid on;