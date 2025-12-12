thesis_startup

num_trials = 100;

for i = 1:2
    if i == 1
        load("../data/processed_data.mat")
        B = data;
        ks = 5:5:80;
    else
        load("../data/large.mat")
        B = full(Problem.A);
        ks = 20:20:400;
    end
    [m,n] = size(B);
    [Z,D,V] = svd(B);
    Bfro = norm(B,"fro");
    
    rpcur2_errs = zeros(length(ks),num_trials);
    rpcurlev_errs = zeros(length(ks),num_trials);
    mdcur_errs = zeros(length(ks),num_trials);
    rpcur2_times = zeros(length(ks),num_trials);
    rpcurlev_times = zeros(length(ks),num_trials);
    mdcur_times = zeros(length(ks),num_trials);
    best = zeros(length(ks));
    
    warning('off', 'all')
    for k_idx = 1:length(ks)
        k = ks(k_idx)
        l = min([ceil(1.5*log(k)*k) m n]);
        kk = ceil(k/(2*log(k)));
        for trial = 1:num_trials
            tic;
            [S,T,~,G,P] = rpcur2(B,k,l);
            rpcur2_times(k_idx,trial) = toc;
            rpcur2_errs(k_idx,trial) =...
                norm(B - B(:,S)*(P\(G*B(T,:))),"fro") / Bfro;

            tic;
            [S,T,~,G,P] = rpcur_lev(B,k,l);
            rpcurlev_times(k_idx,trial) = toc;
            rpcurlev_errs(k_idx,trial) =...
                norm(B - B(:,S)*(P\(G*B(T,:))),"fro") / Bfro;

            if i == 1
                tic
                [S,T,~,G,P] = md_cur(B,k,k,Z(:,1:kk),V(:,1:kk));
                mdcur_times(k_idx,trial) = toc;
            else
                tic
                [ZZ,~,VV] = rsvd(@(X) B*X,@(X) B'*X,size(B,2),kk);
                [S,T,~,G,P] = md_cur(B,k,k,ZZ,VV);
                mdcur_times(k_idx,trial) = toc;
            end
            mdcur_errs(k_idx,trial) =...
                norm(B - B(:,S)*(P\(G*B(T,:))),"fro") / Bfro;
        end
        best(k_idx) = norm(B - Z(:,1:k)*D(1:k,1:k)*V(:,1:k)',"fro") / Bfro;
    end
    warning('on', 'all')
    
    figure
    quantiles = quantile(rpcur2_errs,[0.1 0.5 0.9],2);
    plot_shaded(ks,quantiles(:,2),quantiles(:,1),quantiles(:,3),orange,"Marker","*","MarkerSize",10)
    quantiles = quantile(rpcurlev_errs,[0.1 0.5 0.9],2);
    plot_shaded(ks,quantiles(:,2),quantiles(:,1),quantiles(:,3),blue,"Marker","o","MarkerSize",10,"MarkerFaceColor",blue)
    quantiles = quantile(mdcur_errs,[0.1 0.5 0.9],2);
    plot_shaded(ks,quantiles(:,2),quantiles(:,1),quantiles(:,3),pink,"Marker","x","MarkerSize",10)
    plot(ks,best,"k:")
    set(gca,"YScale","log")
    
    xlabel("Rank $k$")
    ylabel("Relative error $\|\mbox{\boldmath $B$} - \mbox{\boldmath $\widehat{B}$}\|_{\rm F}/\|\mbox{\boldmath $B$}\|_{\rm F}$")
    
    legend({"","RPCUR2",...
            "","RPCURLev",...
            "","MDCUR",...
            "Optimal"}, "Location", "best")

    if i == 1
        exportgraphics(gcf,"../figs/fig_10_4.png")
        saveas(gcf,"../figs/fig_10_4.fig")
    else
        exportgraphics(gcf,"../figs/fig_10_3_a.png")
        saveas(gcf,"../figs/fig_10_3_a.fig")

        figure
        quantiles = quantile(rpcur2_times,[0.1 0.5 0.9],2);
        plot_shaded(ks,quantiles(:,2),quantiles(:,1),quantiles(:,3),orange,"Marker","*","MarkerSize",10)
        quantiles = quantile(rpcurlev_times,[0.1 0.5 0.9],2);
        plot_shaded(ks,quantiles(:,2),quantiles(:,1),quantiles(:,3),blue,"Marker","o","MarkerSize",10,"MarkerFaceColor",blue)
        quantiles = quantile(mdcur_times,[0.1 0.5 0.9],2);
        plot_shaded(ks,quantiles(:,2),quantiles(:,1),quantiles(:,3),pink,"Marker","x","MarkerSize",10)
        set(gca,"YScale","log")
    
        xlabel("Rank $k$")
        ylabel("Runtime (sec)")
    
        drawnow

        exportgraphics(gcf,"../figs/fig_10_3_b.png")
        saveas(gcf,"../figs/fig_10_3_b.fig")
    end
end
