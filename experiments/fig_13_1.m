thesis_startup

n = 1000;
A = rand_with_evals([ones(1,ceil(n/20)) 1e-3*ones(1,n-ceil(n/20))]);
s_list = 20:20:300;
num_trials = 100;

trA = trace(A);
gh_errs = zeros(length(s_list), num_trials);
hutchpp_errs = zeros(length(s_list), num_trials);
hutchpp_res_errs = zeros(length(s_list), num_trials);

for s_idx = 1:length(s_list)
    s = s_list(s_idx);
    for trial = 1:num_trials
        [s trial]
        gh_errs(s_idx,trial) = abs(girard_hutchinson(@(X) A*X,n,s)-trA) / trA;
        hutchpp_errs(s_idx,trial) = abs(hutchpp(@(X) A*X,n,s)-trA) / trA;
        hutchpp_res_errs(s_idx,trial) = abs(hutchpp_resphere(@(X) A*X,n,s)-trA) / trA;
    end
end

%% Plot

quantiles = quantile(gh_errs,[0.1 0.5 0.9],2);
plot_shaded(s_list,quantiles(:,2),quantiles(:,1),quantiles(:,3),yellow,"LineStyle","-.")
quantiles = quantile(hutchpp_errs,[0.1 0.5 0.9],2);
plot_shaded(s_list,quantiles(:,2),quantiles(:,1),quantiles(:,3),blue,"LineStyle","--")
quantiles = quantile(hutchpp_res_errs,[0.1 0.5 0.9],2);
plot_shaded(s_list,quantiles(:,2),quantiles(:,1),quantiles(:,3),purple,"LineStyle","-")
set(gca,"YScale","log")

xlabel("Number of matvecs $s$")
ylabel("Relative error")
legend({"","Girard--Hutchinson","","Hutch++","","Hutch++ (resphered)"},"Location","SouthWest")

exportgraphics(gcf,"../figs/fig_13_1.png")
saveas(gcf,"../figs/fig_13_1.fig")