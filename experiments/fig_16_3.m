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
dA = diag(expm(full(A)));

s_list = 20:20:300;
num_trials = 100;

bks_errs = zeros(length(s_list), num_trials, 2);
udiagpp_errs = zeros(length(s_list), num_trials, 2);
xdiag_errs = zeros(length(s_list), num_trials, 2);
xnysdiag_errs = zeros(length(s_list), num_trials, 2);

figure("Position", [100, 100, 1350, 450])
matvec = @(X) lanczos_fa(@(Z) A*Z,X,40,@exp);

for s_idx = 1:length(s_list)
    s = s_list(s_idx);
    s
    for trial = 1:num_trials
        dest = bks(matvec,n,s);
        bks_errs(s_idx,trial,1)...
            = norm((dest-dA)./dA,"Inf");
        bks_errs(s_idx,trial,2)...
            = norm(dest-dA) / norm(dA);

        dest = udiagpp(matvec,matvec,n,s);
        udiagpp_errs(s_idx,trial,1)...
            = norm((dest-dA)./dA,"Inf");
        udiagpp_errs(s_idx,trial,2)...
            = norm(dest-dA) / norm(dA);

        dest = xdiag(matvec,matvec,n,s);
        xdiag_errs(s_idx,trial,1)...
            = norm((dest-dA) ./ dA,"Inf");
        xdiag_errs(s_idx,trial,2)...
            = norm(dest-dA) / norm(dA);

        dest = xnysdiag(matvec,n,s);
        xnysdiag_errs(s_idx,trial,1)...
            = norm((dest-dA) ./ dA,"Inf");
        xnysdiag_errs(s_idx,trial,2)...
            = norm(dest-dA) / norm(dA);
    end
end

for i = 1:2
    subplot(1,2,i)
    quantiles = quantile(bks_errs(:,:,i),[0.1 0.5 0.9],2);
    plot_shaded(s_list,quantiles(:,2),quantiles(:,1),quantiles(:,3),yellow,"Marker","s","MarkerSize",10)
    quantiles = quantile(udiagpp_errs(:,:,i),[0.1 0.5 0.9],2);
    plot_shaded(s_list,quantiles(:,2),quantiles(:,1),quantiles(:,3),purple,"Marker","*","MarkerSize",10)
    quantiles = quantile(xdiag_errs(:,:,i),[0.1 0.5 0.9],2);
    plot_shaded(s_list,quantiles(:,2),quantiles(:,1),quantiles(:,3),blue,"Marker","x","MarkerSize",10)
    quantiles = quantile(xnysdiag_errs(:,:,i),[0.1 0.5 0.9],2);
    plot_shaded(s_list,quantiles(:,2),quantiles(:,1),quantiles(:,3),orange,"Marker","o","MarkerSize",10)
    set(gca,"YScale","log")
    xlabel("Number of matvecs $s$")
    if i == 1
        ylabel("Maximum relative error")
        legend({"","BKS","","UDiag++","","XDiag","","XNysDiag"},"Location","east")
    else
        ylabel("$\ell_2$ relative error")
    end
end


save("../data/fig_16_3.mat","xnysdiag_errs","xdiag_errs","bks_errs","udiagpp_errs")

exportgraphics(gcf,"../figs/fig_16_3.png")
saveas(gcf,"../figs/fig_16_3.fig")