thesis_startup

n = 1000;
A = rand_with_evals(random_signs(1,n).*(1:n).^(-2));
k = 100;
q_list = 0:2:8;
trials = 100;

rsi_errs = zeros(length(q_list),trials);
rsi_ests = zeros(length(q_list),trials);
rsi_int_errs = zeros(length(q_list),trials);

for q_idx = 1:length(q_list)
    q = q_list(q_idx)
    for trial = 1:trials
        [U,S,V,rsi_ests(q_idx,trial)] = rsi_errest(@(X) A*X, @(X) A*X, n, 100, q);
        rsi_errs(q_idx,trial) = norm(A - U*S*V',"fro");
        [U,S,V] = rsi_int_orth(@(X) A*X, @(X) A*X, n, n, 100, q);
        rsi_int_errs(q_idx,trial) = norm(A - U*S*V',"fro");
    end
end

figure
quantiles = quantile(rsi_ests, [0.1 0.5 0.9], 2);
plot_shaded(q_list,quantiles(:,2),quantiles(:,1),quantiles(:,3),darkblue,"Marker","x","MarkerSize",10,"LineStyle","--")
quantiles = quantile(rsi_errs, [0.1 0.5 0.9], 2);
plot_shaded(q_list,quantiles(:,2),quantiles(:,1),quantiles(:,3),blue,"Marker","x","MarkerSize",10)
quantiles = quantile(rsi_int_errs, [0.1 0.5 0.9], 2);
plot_shaded(q_list,quantiles(:,2),quantiles(:,1),quantiles(:,3),orange,"Marker","o","MarkerSize",10)
set(gca,"YScale","log")

xlabel("Subspace iteration steps $q$")
ylabel("Frobenius-norm error or estimate")

legend({"","Estimate (w/o int. orth.)","","Error (w/o int. orth.)","","Error (w/ int. orth.)"},"Location","northwest")

exportgraphics(gcf,"../figs/fig_20_2.png")
saveas(gcf,"../figs/fig_20_2.fig")
