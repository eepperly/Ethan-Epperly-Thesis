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

bks_errs = zeros(length(s_list), num_trials,2);
xnysdiag_errs = zeros(length(s_list), num_trials,2);
jl_errs = zeros(length(s_list), num_trials,2);
xsymrownorm_errs = zeros(length(s_list), num_trials,2);

matvec_A = @(X) lanczos_fa(@(Z) A*Z,X,40,@exp);
matvec_Ahalf = @(X) lanczos_fa(@(Z) A*Z,X,40,@(d) exp(d/2));

for s_idx = 1:length(s_list)
    s = s_list(s_idx);
    s
    for trial = 1:num_trials
        dest = bks(matvec_A,n,s);
        bks_errs(s_idx,trial,1)...
            = norm((dest-dA)./dA,"Inf");
        bks_errs(s_idx,trial,2)...
            = norm(dest-dA) / norm(dA);

        dest = xnysdiag(matvec_A,n,s);
        xnysdiag_errs(s_idx,trial,1)...
            = norm((dest-dA)./dA,"Inf");
        xnysdiag_errs(s_idx,trial,2)...
            = norm(dest-dA) / norm(dA);

        dest = jl_rownorm(matvec_Ahalf,n,s);
        jl_errs(s_idx,trial,1)...
            = norm((dest-dA) ./ dA,"Inf");
        jl_errs(s_idx,trial,2)...
            = norm(dest-dA) / norm(dA);

        dest = xsymrownorm(matvec_Ahalf,n,s);
        xsymrownorm_errs(s_idx,trial,1)...
            = norm((dest-dA) ./ dA,"Inf");
        xsymrownorm_errs(s_idx,trial,2)...
            = norm(dest-dA) / norm(dA);
    end
end

%% Plot

figure("Position", [100, 100, 1350, 450])

for i = 1:2
    subplot(1,2,i)
    quantiles = quantile(bks_errs(:,:,i),[0.1 0.5 0.9],2);
    plot_shaded(s_list,quantiles(:,2),quantiles(:,1),quantiles(:,3),yellow,"Marker","s","MarkerSize",10)
    quantiles = quantile(xnysdiag_errs(:,:,i),[0.1 0.5 0.9],2);
    plot_shaded(s_list,quantiles(:,2),quantiles(:,1),quantiles(:,3),purple,"Marker","*","MarkerSize",10)
    quantiles = quantile(jl_errs(:,:,i),[0.1 0.5 0.9],2);
    plot_shaded(s_list,quantiles(:,2),quantiles(:,1),quantiles(:,3),blue,"Marker","x","MarkerSize",10)
    quantiles = quantile(xsymrownorm_errs(:,:,i),[0.1 0.5 0.9],2);
    plot_shaded(s_list,quantiles(:,2),quantiles(:,1),quantiles(:,3),orange,"Marker","o","MarkerSize",10)
    set(gca,"YScale","log")
    xlabel("Number of matvecs $s$")
    if i == 1
        title("max error","FontName","Courier New","Interpreter","LaTeX")
        ylabel("Maximum relative error")
        legend({"","BKS","","XNysDiag","","JL","","XSymRowNorm"},"Location","east")
    else
        title("$\ell_2$ error","FontName","Courier New","Interpreter","LaTeX")
        ylabel("$\ell_2$ relative error")
    end
end

save("../data/fig_17_3.mat","xnysdiag_errs","bks_errs","jl_errs","xsymrownorm_errs")

exportgraphics(gcf,"../figs/fig_17_3.png")
saveas(gcf,"../figs/fig_17_3.fig")