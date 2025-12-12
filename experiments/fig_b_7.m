ns = round(logspace(1,log10(600),11));
ms = ns.^2;

trials = 100;

sparse_times = zeros(length(ns), trials, 4);
gaussian_times = zeros(length(ns), trials, 4);
srtt_times = nan*zeros(length(ns), trials, 4);

for n_idx = 1:length(ns)
    m = ms(n_idx); n = ns(n_idx); d = 2*n;
    n
    B = randn(m,n);
    c = randn(m,1);
    F = sparse_sign(n,m,10)';
    for trial = 1:trials
        tic; St = randn(d,m)/sqrt(d); gaussian_times(n_idx,trial,1) = toc;
        tic; St*c; gaussian_times(n_idx,trial,2) = toc;
        tic; St*B; gaussian_times(n_idx,trial,3) = toc;
        tic; St*F; gaussian_times(n_idx,trial,4) = toc;

        tic; St = sparse_sign(d,m,4); sparse_times(n_idx,trial,1) = toc;
        tic; St*c; sparse_times(n_idx,trial,2) = toc;
        tic; St*B; sparse_times(n_idx,trial,3) = toc;
        tic; St*F; sparse_times(n_idx,trial,4) = toc;

        tic; St = srtt(d,m,1); srtt_times(n_idx,trial,1) = toc;
        tic; St*c; srtt_times(n_idx,trial,2) = toc;
        tic; St*B; srtt_times(n_idx,trial,3) = toc;
    end
end

%% Plot

figure("Position", [100, 100, 1350, 900])

for i = 1:4
    subplot(2,2,i)
    
    quantiles = quantile(gaussian_times(:,:,i),[0.1 0.5 0.9],2);
    plot_shaded(ns,quantiles(:,2),quantiles(:,1),quantiles(:,3),purple,"Marker","*","MarkerSize",10)
    quantiles = quantile(srtt_times(:,:,i),[0.1 0.5 0.9],2);
    plot_shaded(ns,quantiles(:,2),quantiles(:,1),quantiles(:,3),blue,"Marker","x","MarkerSize",10)
    quantiles = quantile(sparse_times(:,:,i),[0.1 0.5 0.9],2);
    plot_shaded(ns,quantiles(:,2),quantiles(:,1),quantiles(:,3),orange,"Marker","o","MarkerSize",10)
    set(gca,"YScale","log")
    set(gca,"XScale","log")

    xlabel("Dimension $n = m^{1/2}$")
    ylabel("Time (sec)")

    axis([-Inf Inf 1e-6 1e2])

    if i == 1
        legend({"","Gaussian","","SRTT","","Sparse"},"Location","northwest")
    end
    drawnow
end

subplot(2,2,1)
title("initialize","FontName","Courier New","Interpreter","TeX")
subplot(2,2,2)
title("vector apply","FontName","Courier New","Interpreter","TeX")
subplot(2,2,3)
title("dense apply","FontName","Courier New","Interpreter","TeX")
subplot(2,2,4)
title("sparse apply","FontName","Courier New","Interpreter","TeX")

exportgraphics(gcf,"../figs/fig_b_7.png")
saveas(gcf,"../figs/fig_b_7.fig")
