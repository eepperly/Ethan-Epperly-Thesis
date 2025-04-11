thesis_startup

n = 1000;
As = {rand_with_evals(linspace(1,3,n)),...
      rand_with_evals((1:n).^(-2)),...
      rand_with_evals((0.7.^(0:(n-1)))),...
      rand_with_evals([ones(1,ceil(n/20)) 1e-3*ones(1,n-ceil(n/20))])};
names = {'flat','poly','exp','step'}; 

s_list = 20:20:300;
num_trials = 1000;

gh_errs = zeros(length(s_list), num_trials, 4);
hutchpp_errs = zeros(length(s_list), num_trials, 4);
xtrace_errs = zeros(length(s_list), num_trials, 4);
xnystrace_errs = zeros(length(s_list), num_trials, 4);

figure("Position", [100, 100, 1350, 900])

for A_idx = 1:length(As)
    A = As{A_idx};
    trA = trace(A);

    for s_idx = 1:length(s_list)
        s = s_list(s_idx);
        [A_idx s]
        for trial = 1:num_trials
            gh_errs(s_idx,trial,A_idx)...
                = abs(girard_hutchinson(@(X) A*X,n,s)-trA) / abs(trA);
            hutchpp_errs(s_idx,trial,A_idx)...
                = abs(hutchpp(@(X) A*X,n,s)-trA) / abs(trA);
            xtrace_errs(s_idx,trial,A_idx)...
                = abs(xtrace(@(X) A*X,n,s)-trA) / abs(trA);
            xnystrace_errs(s_idx,trial,A_idx)...
                = abs(xnystrace(@(X) A*X,n,s)-trA) / abs(trA);
        end
    end

    subplot(2,2,A_idx)
    quantiles = quantile(gh_errs(:,:,A_idx),[0.1 0.5 0.9],2);
    plot_shaded(s_list,quantiles(:,2),quantiles(:,1),quantiles(:,3),yellow,"Marker","s","MarkerSize",10)
    quantiles = quantile(hutchpp_errs(:,:,A_idx),[0.1 0.5 0.9],2);
    plot_shaded(s_list,quantiles(:,2),quantiles(:,1),quantiles(:,3),purple,"Marker","*","MarkerSize",10)
    quantiles = quantile(xtrace_errs(:,:,A_idx),[0.1 0.5 0.9],2);
    plot_shaded(s_list,quantiles(:,2),quantiles(:,1),quantiles(:,3),blue,"Marker","x","MarkerSize",10)
    quantiles = quantile(xnystrace_errs(:,:,A_idx),[0.1 0.5 0.9],2);
    plot_shaded(s_list,quantiles(:,2),quantiles(:,1),quantiles(:,3),orange,"Marker","o","MarkerSize",10)
    set(gca,"YScale","log")

    xlabel("Number of matvecs $s$")
    ylabel("Relative error")

    if A_idx == 3
        legend({"","Girard--Hutchinson","","Hutch++","","XTrace","","XNysTrace"},"Location","east")
    end
    drawnow
end

save("../data/fig_15_1.mat","xnystrace_errs","xtrace_errs","gh_errs","hutchpp_errs")

subplot(2,2,1)
title("flat","FontName","Courier New","Interpreter","TeX")
subplot(2,2,2)
title("poly","FontName","Courier New","Interpreter","TeX")
subplot(2,2,3)
title("exp","FontName","Courier New","Interpreter","TeX")
subplot(2,2,4)
title("step","FontName","Courier New","Interpreter","TeX")

exportgraphics(gcf,"../figs/fig_15_1.png")
saveas(gcf,"../figs/fig_15_1.fig")