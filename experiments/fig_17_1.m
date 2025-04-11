thesis_startup

n = 1000;
As = {rand_with_evals(linspace(1,3,n)),...
      rand_with_evals((1:n).^(-2)),...
      rand_with_evals((0.7.^(0:(n-1)))),...
      rand_with_evals([ones(1,ceil(n/20)) 1e-3*ones(1,n-ceil(n/20))])};
names = {'flat','poly','exp','step'}; 

s_list = 20:20:300;
num_trials = 100;

bks_errs = zeros(length(s_list), num_trials, 4);
diagpp_errs = zeros(length(s_list), num_trials, 4);
xdiag_errs = zeros(length(s_list), num_trials, 4);
xnysdiag_errs = zeros(length(s_list), num_trials, 4);

figure("Position", [100, 100, 1350, 900])

for A_idx = 1:length(As)
    A = As{A_idx};
    dA = diag(A);

    for s_idx = 1:length(s_list)
        s = s_list(s_idx);
        [A_idx s]
        for trial = 1:num_trials
            bks_errs(s_idx,trial,A_idx)...
                = norm(bks(@(X) A*X,n,s)-dA,Inf) / norm(dA,Inf);
            diagpp_errs(s_idx,trial,A_idx)...
                = norm(udiagpp(@(X) A*X,@(X) A*X,n,s)-dA,Inf) / norm(dA,Inf);
            xdiag_errs(s_idx,trial,A_idx)...
                = norm(xdiag(@(X) A*X,@(X) A*X,n,s)-dA,Inf) / norm(dA,Inf);
            xnysdiag_errs(s_idx,trial,A_idx)...
                = norm(xnysdiag(@(X) A*X,n,s)-dA,Inf) / norm(dA,Inf);
        end
    end

    subplot(2,2,A_idx)
    quantiles = quantile(bks_errs(:,:,A_idx),[0.1 0.5 0.9],2);
    plot_shaded(s_list,quantiles(:,2),quantiles(:,1),quantiles(:,3),yellow,"Marker","s","MarkerSize",10)
    quantiles = quantile(diagpp_errs(:,:,A_idx),[0.1 0.5 0.9],2);
    plot_shaded(s_list,quantiles(:,2),quantiles(:,1),quantiles(:,3),purple,"Marker","*","MarkerSize",10)
    quantiles = quantile(xdiag_errs(:,:,A_idx),[0.1 0.5 0.9],2);
    plot_shaded(s_list,quantiles(:,2),quantiles(:,1),quantiles(:,3),blue,"Marker","x","MarkerSize",10)
    quantiles = quantile(xnysdiag_errs(:,:,A_idx),[0.1 0.5 0.9],2);
    plot_shaded(s_list,quantiles(:,2),quantiles(:,1),quantiles(:,3),orange,"Marker","o","MarkerSize",10)
    set(gca,"YScale","log")

    xlabel("Number of matvecs $s$")
    ylabel("Relative error")

    if A_idx == 3
        legend({"","BKS","","UDiag++","","XDiag","","XNysDiag"},"Location","east")
    end
    drawnow
end

save("../data/fig_15_1.mat","xnysdiag_errs","xdiag_errs","bks_errs","diagpp_errs")

subplot(2,2,1)
title("flat","FontName","Courier New","Interpreter","TeX")
subplot(2,2,2)
title("poly","FontName","Courier New","Interpreter","TeX")
subplot(2,2,3)
title("exp","FontName","Courier New","Interpreter","TeX")
subplot(2,2,4)
title("step","FontName","Courier New","Interpreter","TeX")

exportgraphics(gcf,"../figs/fig_17_1.png")
saveas(gcf,"../figs/fig_17_1.fig")