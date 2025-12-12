thesis_startup

m = 1e5;
n = 1e3;
xx = linspace(-1,1,m).';
B = cos(acos(xx) * (0:(n-1)));
Q = orth(B);

srnQ = sqrownorms(Q);

s_list = 4:4:60;
num_trials = 100;

jl_errs = zeros(length(s_list), num_trials);
bks_errs = zeros(length(s_list), num_trials);
xnysdiag_errs = zeros(length(s_list), num_trials);

for s_idx = 1:length(s_list)
    s = s_list(s_idx)
    for trial = 1:num_trials
        jl_est = jl_rownorm(@(X) Q*X,n,s);
        jl_errs(s_idx,trial) = max([jl_est./srnQ;srnQ./jl_est]);
        bks_est = (bks(@(X) Q*(Q'*X),m,s/2)-srnQ)./srnQ;
        bks_errs(s_idx,trial) = max([bks_est./srnQ;srnQ./bks_est]);
        xnysdiag_est = xnysdiag(@(X) Q*(Q'*X),m,s/2);
        xnysdiag_errs(s_idx,trial) = max([xnysdiag_est./srnQ;srnQ./xnysdiag_est]);
    end
end

%% Plot

figure(1)

quantiles = quantile(xnysdiag_errs,[0.1 0.5 0.9],2);
plot_shaded(s_list,quantiles(:,2),quantiles(:,1),quantiles(:,3),purple,"Marker","*","MarkerSize",10)
quantiles = quantile(bks_errs,[0.1 0.5 0.9],2);
plot_shaded(s_list,quantiles(:,2),quantiles(:,1),quantiles(:,3),yellow,"Marker","s","MarkerSize",10)
quantiles = quantile(jl_errs,[0.1 0.5 0.9],2);
plot_shaded(s_list,quantiles(:,2),quantiles(:,1),quantiles(:,3),orange,"Marker","o","MarkerSize",10)

set(gca,"YScale","log")

xlabel("Number of matvecs $s$")
ylabel("Approximation factor")

legend({"","XNysDiag","","BKS","","JL"},"Location","east")

exportgraphics(gcf,"../figs/fig_17_1_a.png")
saveas(gcf,"../figs/fig_17_1_a.fig")

figure(2)

alphaVal = 0.1;
scatter(xx, xnysdiag_est, 100, hex2rgb(purple), '*', ...
    'MarkerFaceAlpha', alphaVal, 'MarkerEdgeAlpha', alphaVal); hold on
scatter(xx, bks_est, 100, hex2rgb(yellow), 's', ...
    'MarkerFaceAlpha', alphaVal, 'MarkerEdgeAlpha', alphaVal);
scatter(xx, jl_est, 100, hex2rgb(orange), 'o', ...
    'MarkerFaceAlpha', alphaVal, 'MarkerEdgeAlpha', alphaVal);semilogy(xx,srnQ,"Color",black,"LineStyle","--","LineWidth",4);
plot(xx, srnQ, "Color", black, "LineStyle", "--", "LineWidth", 4);

set(gca,"YScale","log")

ylabel("Leverage score")

axis([-1 1 1e-4 1e1])

exportgraphics(gcf,"../figs/fig_17_1_b.png")
saveas(gcf,"../figs/fig_17_1_b.fig")