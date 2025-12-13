thesis_startup

n = 1000;
As = {rand_with_evals(random_signs(1,n).*linspace(1,3,n)),...
      rand_with_evals(random_signs(1,n).*(1:n).^(-2)),...
      rand_with_evals(random_signs(1,n).*(0.7.^(0:(n-1)))),...
      rand_with_evals(random_signs(1,n).*...
      [ones(1,ceil(n/20)) 1e-3*ones(1,n-ceil(n/20))])};
names = {'flat','poly','exp','step'}; 

s_list = 20:20:300;
num_trials = 100;

figure("Position", [100, 100, 900, 600])

for A_idx = 1:length(As)
    A = As{A_idx};
    trA = trace(A);

    xtrace_errs = zeros(length(s_list), num_trials);
    xsymtrace_errs = zeros(length(s_list), num_trials);

    for s_idx = 1:length(s_list)
        s = s_list(s_idx);
        for trial = 1:num_trials
            xtrace_errs(s_idx,trial)...
                = abs(xtrace(@(X) A*X,n,s)-trA) / abs(trA);
            xsymtrace_errs(s_idx,trial)...
                = abs(xsymtrace(@(X) A*X,n,s)-trA) / abs(trA);
        end
    end

    subplot(2,2,A_idx)
    quantiles = quantile(xtrace_errs,[0.2 0.5 0.8],2);
    plot_shaded(s_list,quantiles(:,2),quantiles(:,1),quantiles(:,3),orange)
    quantiles = quantile(xsymtrace_errs,[0.2 0.5 0.8],2);
    plot_shaded(s_list,quantiles(:,2),quantiles(:,1),quantiles(:,3),blue,"LineStyle","--")
    set(gca,"YScale","log")

    xlabel("Number of matvecs $s$")
    ylabel("Relative error")

    if A_idx == 3
        legend("","XTrace","","XSymTrace")
    end
end

subplot(2,2,1)
title("flat","FontName","Courier New","Interpreter","TeX")
subplot(2,2,2)
title("poly","FontName","Courier New","Interpreter","TeX")
subplot(2,2,3)
title("exp","FontName","Courier New","Interpreter","TeX")
subplot(2,2,4)
title("step","FontName","Courier New","Interpreter","TeX")

exportgraphics(gcf,"../figs/fig_20_1.png")
saveas(gcf,"../figs/fig_20_1.fig")