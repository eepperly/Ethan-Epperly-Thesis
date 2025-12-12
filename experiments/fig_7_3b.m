thesis_startup

trials = 100;
ks = 5:5:60;
kernel = @(S,T) exp(-pdist2(S,T).^2/2);
f = @(S) cos(sqrownorms(S).^0.5);
intf = 0.2752215409929235;
% f = @(S) ones(size(S,1),1);
% intf = 1;

rpc_errs = zeros(length(ks), trials);
unif_errs = zeros(length(ks), trials);
mc_errs = zeros(length(ks), trials);

for i = 1:length(ks)
    k = ks(i)
    for trial = 1:trials
        S = rpcholesky_rejection(kernel,@() randn(1,2),k,2);
        Au = exp(-sqrownorms(S)/4) / 2;
        [~,integrator] = kernel_quad_wts(kernel,S,Au);
        rpc_errs(i,trial) = abs(integrator(f) - intf) / intf;

        S = randn(k,2);
        Au = exp(-sqrownorms(S)/4) / 2;
        [~,integrator] = kernel_quad_wts(kernel,S,Au);
        unif_errs(i,trial) = abs(integrator(f) - intf) / intf;
        mc_errs(i,trial) = abs(mean(f(S)) - intf) / intf;
    end
end

unif_means = mean(unif_errs,2);
unif_stds = std(unif_errs,0,2);
rpc_means = mean(rpc_errs,2);
rpc_stds = std(rpc_errs,0,2);
mc_means = mean(mc_errs,2);
mc_stds = std(mc_errs,0,2);

figure
plot(ks, mc_means, "Color", purple, "LineStyle", "-."); hold on
plot_shaded(ks, unif_means, unif_means-unif_stds, unif_means+unif_stds, blue, "LineStyle", "--"); hold on
plot_shaded(ks, rpc_means, rpc_means-rpc_stds, rpc_means+rpc_stds, orange)
set(gca,"YScale","log")

legend({"Monte Carlo","","IID","","RPCholesky"},'FontSize',25,'Location','Southwest')
xlabel("Number of landmarks $k$")
ylabel("Relative quadrature error")

saveas(gcf,"../figs/fig_7_3_b.png")
saveas(gcf,"../figs/fig_7_3_b.fig")

function s = propose_special_region()
while true
    r = 2*rand()^(1/4);
    theta = 2*pi*rand();
    x = r*cos(theta);
    y = r*sin(theta);
    if (x-1)^2 + y^2 > 1
        break
    end
end
s = [x;y];
end

function K = matern32(S,T)
    X = sqrt(3)*pdist2(S,T);
    K = (1 + X) .* exp(-X);
end