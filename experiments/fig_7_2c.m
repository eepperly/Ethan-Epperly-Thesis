thesis_startup

f = @(S) 1 ./ (1 + sqrownorms(S).^0.5); % f(x) = 1/(1+||x||)
laplace = @(S,T) exp(-pdist2(S,T));

ks = 10:10:300;
mctrials = 1e3;
trials = 100;
unif_mses = zeros(length(ks),trials);
rpc_mses = zeros(length(ks),trials);
for i = 1:length(ks)
    k = ks(i)
    for trial = 1:trials
        % Uniform 
        S = randn(k,2);
        fhat = kernel_interp(f,laplace,S);
        U = randn(mctrials,2);
        unif_mses(i,trial) = mean((f(U) - fhat(U)).^2);
    
        % RPCholesky 
        S = rpcholesky_rejection(laplace,@() randn(1,2),k,2);
        fhat = kernel_interp(f,laplace,S);
        U = randn(mctrials,2);
        rpc_mses(i,trial) = mean((f(U) - fhat(U)).^2);
    end
end

ts = 3:17;
ten_mses = zeros(length(ts),trials);
for i = 1:length(ts)
    t = ts(i)
    for trial = 1:trials
        X = norminv(linspace(1/(2*t),1-1/(2*t),t));
        [xx,yy] = meshgrid(X);
        S = [xx(:) yy(:)];
        fhat = kernel_interp(f,laplace,S);
        U = randn(mctrials,2);
        ten_mses(i,trial) = mean((f(U) - fhat(U)).^2);
    end
end

unif_means = mean(unif_mses,2);
unif_stds = std(unif_mses,0,2);
rpc_means = mean(rpc_mses,2);
rpc_stds = std(rpc_mses,0,2);
ten_means = mean(ten_mses,2);
ten_stds = std(ten_mses,0,2);

%% Plot

figure('Position', [100, 100, 700, 400]);

plot_shaded(ts.^2, ten_means, ten_means-ten_stds, ten_means+ten_stds, purple, "LineStyle", "-.")
plot_shaded(ks, unif_means, unif_means-unif_stds, unif_means+unif_stds, blue, "LineStyle", "--"); hold on
plot_shaded(ks, rpc_means, rpc_means-rpc_stds, rpc_means+rpc_stds, orange)
set(gca,"YScale","log")

legend("","Tensor","","Uniform","","RPCholesky")
xlabel("Number of landmarks $k$")
ylabel("Mean-square error")

saveas(gcf,"../figs/fig_7_2_c.png")
saveas(gcf,"../figs/fig_7_2_c.fig")