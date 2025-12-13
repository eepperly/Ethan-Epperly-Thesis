thesis_startup

try 
    load("../data/yeast.mat")
catch
    cd("../data")
    ssget("Pajek/yeast")
    cd("../experiments")
end

A = Problem.A;
n = size(A,1);
tr = trace(expm(full(A)));

s_list = 20:20:300;
num_trials = 100;

gh_errs = zeros(length(s_list), num_trials);
hutchpp_errs = zeros(length(s_list), num_trials);
xtrace_errs = zeros(length(s_list), num_trials);
xnystrace_errs = zeros(length(s_list), num_trials);

figure

for s_idx = 1:length(s_list)
    s = s_list(s_idx);
    s
    for trial = 1:num_trials
        gh_errs(s_idx,trial)...
            = abs(girard_hutchinson(@(X) lanczos_fa(@(Z) A*Z,...
            X,40,@exp),n,s)-tr) / tr;
        hutchpp_errs(s_idx,trial)...
            = abs(hutchpp_resphere(@(X) lanczos_fa(@(Z) A*Z,...
            X,40,@exp),n,s)-tr) / tr;
        xtrace_errs(s_idx,trial)...
            = abs(xtrace_resphere(@(X) lanczos_fa(@(Z) A*Z,...
            X,40,@exp),n,s)-tr) / tr;
        xnystrace_errs(s_idx,trial)...
            = abs(xnystrace_resphere(@(X) lanczos_fa(@(Z) A*Z,...
            X,40,@exp),n,s)-tr) / tr;
    end
end

quantiles = quantile(gh_errs,[0.1 0.5 0.9],2);
plot_shaded(s_list,quantiles(:,2),quantiles(:,1),quantiles(:,3),yellow,"Marker","s","MarkerSize",10)
quantiles = quantile(hutchpp_errs,[0.1 0.5 0.9],2);
plot_shaded(s_list,quantiles(:,2),quantiles(:,1),quantiles(:,3),purple,"Marker","*","MarkerSize",10)
quantiles = quantile(xtrace_errs,[0.1 0.5 0.9],2);
plot_shaded(s_list,quantiles(:,2),quantiles(:,1),quantiles(:,3),blue,"Marker","x","MarkerSize",10)
quantiles = quantile(xnystrace_errs,[0.1 0.5 0.9],2);
plot_shaded(s_list,quantiles(:,2),quantiles(:,1),quantiles(:,3),orange,"Marker","o","MarkerSize",10)
set(gca,"YScale","log")

xlabel("Number of matvecs $s$")
ylabel("Relative error")

legend({"","Girard--Hutchinson","","Hutch++","","XTrace","","XNysTrace"},"Location","east")

save("../data/fig_14_4.mat","xnystrace_errs","xtrace_errs","gh_errs","hutchpp_errs")


exportgraphics(gcf,"../figs/fig_14_4.png")
saveas(gcf,"../figs/fig_14_4.fig")
