thesis_startup

n = 1000;
As = {rand_with_evals(linspace(1,3,n)),...
      rand_with_evals((1:n).^(-2)),...
      rand_with_evals((0.7.^(0:(n-1)))),...
      rand_with_evals([ones(1,ceil(n/20)) 1e-3*ones(1,n-ceil(n/20))])};
names = {'flat','poly','exp','step'}; 

s_list = 20:20:300;
num_trials = 1000;

xtrace_errs = zeros(length(s_list), num_trials, 4);
xtrace_ests = zeros(length(s_list), num_trials, 4);

figure("Position", [100, 100, 1350, 900])

for A_idx = 1:length(As)
    A = As{A_idx};
    trA = trace(A);

    for s_idx = 1:length(s_list)
        s = s_list(s_idx)
        for trial = 1:num_trials
            [est,err] = xtrace(@(X) A*X,n,s);
            rel_err = abs(est - trA) / abs(trA);
            xtrace_errs(s_idx,trial,A_idx) = rel_err;
            xtrace_ests(s_idx,trial,A_idx) = err / trA;
        end
    end

    subplot(2,2,A_idx)
    
    % Plot on the left y-axis: relative error and its estimate
    yyaxis left
    plot(s_list,median(xtrace_errs(:,:,A_idx),2),"Color","#0072BD","Marker","x","MarkerSize",10); hold on
    plot(s_list,median(xtrace_ests(:,:,A_idx),2),"--","Color","#0072BD","Marker","x","MarkerSize",10)
    plot(nan,":","Color","#D95319","LineWidth",4)
    ylabel("Relative error or estimate")
    set(gca,"YScale","log")

    % Plot on the right y-axis: error-to-estimate ratio
    yyaxis right
    plot(s_list,median(xtrace_ests(:,:,A_idx),2) ./ median(xtrace_errs(:,:,A_idx),2),":","Color","#D95319","LineWidth",4)
    ylabel("Estimate-to-error ratio")

    xlabel("Number of matvecs $s$")

    if A_idx == 3
        yyaxis left
        legend({"Error","Estimate","Ratio"},"Location","northeast")
    end

    drawnow
end

subplot(2,2,1)
title("flat","FontName","Courier New","Interpreter","TeX")
subplot(2,2,2)
title("poly","FontName","Courier New","Interpreter","TeX")
subplot(2,2,3)
title("exp","FontName","Courier New","Interpreter","TeX")
subplot(2,2,4)
title("step","FontName","Courier New","Interpreter","TeX")

exportgraphics(gcf,"../figs/fig_15_2.png")
saveas(gcf,"../figs/fig_15_2.fig")