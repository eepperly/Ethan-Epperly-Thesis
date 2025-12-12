thesis_startup

n = 1000;
As = {rand_with_evals(linspace(1,3,n)),...
      rand_with_evals((1:n).^(-2)),...
      rand_with_evals((0.7.^(0:(n-1)))),...
      rand_with_evals([ones(1,ceil(n/20)) 1e-3*ones(1,n-ceil(n/20))])};
names = {'flat','poly','exp','step'}; 

s_list = 20:20:300;
num_trials = 100;

jl_errs = zeros(length(s_list), num_trials, 4);
sl3_errs = zeros(length(s_list), num_trials, 4);
sl4_errs = zeros(length(s_list), num_trials, 4);
xrownorm_errs = zeros(length(s_list), num_trials, 4);
xsymrownorm_errs = zeros(length(s_list), num_trials, 4);

figure("Position", [100, 100, 1350, 900])

for A_idx = 1:length(As)
    A = As{A_idx};
    srnA = sqrownorms(A);

    for s_idx = 1:length(s_list)
        s = s_list(s_idx);
        [A_idx s]
        for trial = 1:num_trials
            jl_errs(s_idx,trial,A_idx)...
                = norm((jl_rownorm(@(X) A*X,n,s)-srnA)./srnA,Inf);
            sl4_errs(s_idx,trial,A_idx)...
                = norm((sl4(@(X) A*X,@(X) A*X,n,n,s)-srnA)./srnA,Inf);
            sl3_errs(s_idx,trial,A_idx)...
                = norm((sl3(@(X) A*X,@(X) A*X,n,n,s)-srnA)./srnA,Inf);
            xrownorm_errs(s_idx,trial,A_idx)...
                = norm((xrownorm(@(X) A*X,@(X) A*X,n,s)-srnA)./srnA,Inf);
            xsymrownorm_errs(s_idx,trial,A_idx)...
                = norm((xsymrownorm(@(X) A*X,n,s)-srnA)./srnA,Inf);
        end
    end

    subplot(2,2,A_idx)
    quantiles = quantile(jl_errs(:,:,A_idx),[0.1 0.5 0.9],2);
    plot_shaded(s_list,quantiles(:,2),quantiles(:,1),quantiles(:,3),yellow,"Marker","s","MarkerSize",10)
    quantiles = quantile(sl4_errs(:,:,A_idx),[0.1 0.5 0.9],2);
    plot_shaded(s_list,quantiles(:,2),quantiles(:,1),quantiles(:,3),purple,"Marker","*","MarkerSize",10)
    quantiles = quantile(sl3_errs(:,:,A_idx),[0.1 0.5 0.9],2);
    plot_shaded(s_list,quantiles(:,2),quantiles(:,1),quantiles(:,3),pink,"Marker","^","MarkerSize",10)
    quantiles = quantile(xrownorm_errs(:,:,A_idx),[0.1 0.5 0.9],2);
    plot_shaded(s_list,quantiles(:,2),quantiles(:,1),quantiles(:,3),blue,"Marker","x","MarkerSize",10)
    quantiles = quantile(xsymrownorm_errs(:,:,A_idx),[0.1 0.5 0.9],2);
    plot_shaded(s_list,quantiles(:,2),quantiles(:,1),quantiles(:,3),orange,"Marker","o","MarkerSize",10)
    set(gca,"YScale","log")

    xlabel("Number of matvecs $s$")
    ylabel("Maximum relative error")

    if A_idx == 3
        legend({"","JL","","SL4","","SL3","","XRowNorm","","XSymRowNorm"},"Location","east")
    end
    drawnow
end

save("../data/fig_17_2.mat","xsymrownorm_errs","xrownorm_errs","jl_errs","sl3_errs","sl4_errs")

subplot(2,2,1)
title("flat","FontName","Courier New","Interpreter","TeX")
subplot(2,2,2)
title("poly","FontName","Courier New","Interpreter","TeX")
subplot(2,2,3)
title("exp","FontName","Courier New","Interpreter","TeX")
subplot(2,2,4)
title("step","FontName","Courier New","Interpreter","TeX")

exportgraphics(gcf,"../figs/fig_17_2.png")
saveas(gcf,"../figs/fig_17_2.fig")