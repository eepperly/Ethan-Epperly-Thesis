thesis_startup

n = 1000;
A = rand_with_evals((0.7.^(0:(n-1))));

s_list = 20:20:300;
num_trials = 1000;

xnystrace_bad_errs = zeros(length(s_list), 4);

figure
trA = trace(A);

for s_idx = 1:length(s_list)
    s = s_list(s_idx);
    s
    for trial = 1:num_trials
        xnystrace_bad_errs(s_idx,trial)...
            = abs(xnystrace(@(X) A*X,n,s)-trA) / abs(trA);
    end
end

load("../data/fig_15_1.mat")

%% Plot

quantiles = quantile(xnystrace_bad_errs,[0.1 0.5 0.9],2);
plot_shaded(s_list,quantiles(:,2),quantiles(:,1),quantiles(:,3),pink,"Marker","^","MarkerSize",10,"LineStyle","--"); hold on
quantiles = quantile(xnystrace_errs(:,:,3),[0.1 0.5 0.9],2);
plot_shaded(s_list,quantiles(:,2),quantiles(:,1),quantiles(:,3),orange,"Marker","o","MarkerSize",10)

set(gca,"YScale","log")

xlabel("Number of matvecs $s$")
ylabel("Relative error")

legend({"","MacOS/Apple Silicon","","Linux/x86"},"Location","NorthEast")

exportgraphics(gcf,"../figs/fig_15_2.png")
saveas(gcf,"../figs/fig_15_2.fig")