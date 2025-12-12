thesis_startup

n = 1000;
As = {rand_with_evals(linspace(1,3,n)),...
      rand_with_evals((1:n).^(-2)),...
      rand_with_evals((0.7.^(0:(n-1)))),...
      rand_with_evals([ones(1,ceil(n/20)) 1e-3*ones(1,n-ceil(n/20))])};
names = {'flat','poly','exp','step'}; 

k_list = 10:10:150;
num_trials = 1000;

leaveout_errs = zeros(length(k_list), num_trials, 4);
gh_errs = zeros(length(k_list), num_trials, 4);

figure("Position", [100, 100, 1350, 900])

for A_idx = 1:length(As)
    A = As{A_idx};

    for k_idx = 1:length(k_list)
        k = k_list(k_idx);
        [A_idx k]
        for trial = 1:num_trials
            [U,S,V,est] = rsvd_errest(@(X) A*X,@(X) A*X,n,k);
            err = norm(A-U*S*V',"fro");
            leaveout_errs(k_idx,trial) = max(est / err,err / est) - 1;
            Om = random_signs(1000,10);
            est = norm(A*Om - U*(S*(V'*Om)),"fro")/sqrt(10);
            gh_errs(k_idx,trial) = max(est / err,err / est) - 1;
        end
    end

    subplot(2,2,A_idx)
    quantiles = quantile(gh_errs,[0.1 0.5 0.9],2);
    plot_shaded(k_list,quantiles(:,2),quantiles(:,1),quantiles(:,3),purple,"Marker","*","MarkerSize",10)
    quantiles = quantile(leaveout_errs,[0.1 0.5 0.9],2);
    plot_shaded(k_list,quantiles(:,2),quantiles(:,1),quantiles(:,3),orange,"Marker","o","MarkerSize",10)
    set(gca,"YScale","log")

    xlabel("Approximation rank $k$")
    ylabel("Approximation factor")

    if A_idx == 2
        legend({"","Girard--Hutchinson","","Leave-one-out"},"Location","southwest")
    end
    drawnow
end

save("../data/fig_18_1.mat","gh_errs","leaveout_errs")

subplot(2,2,1)
title("flat","FontName","Courier New","Interpreter","TeX")
subplot(2,2,2)
title("poly","FontName","Courier New","Interpreter","TeX")
subplot(2,2,3)
title("exp","FontName","Courier New","Interpreter","TeX")
subplot(2,2,4)
title("step","FontName","Courier New","Interpreter","TeX")

exportgraphics(gcf,"../figs/fig_18_1.png")
saveas(gcf,"../figs/fig_18_1.fig")