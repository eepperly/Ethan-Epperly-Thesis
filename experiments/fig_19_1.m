thesis_startup

n = 1000;
A = rand_with_evals((0.7.^(0:(n-1))));

trials = 100;
k_list = 25:25:150;

fs = {@(d) [ones(10,1);zeros(length(d)-10,1)], @sqrt};

jacks = zeros(length(k_list),trials,2);
stds = zeros(length(k_list),1,2);
errs = zeros(length(k_list),trials,2);

figure("Position", [100, 100, 1350, 450])

for f_idx = 1:2
    f = fs{f_idx};
    [Q,d] = eig(A,"vector");
    d = d(end:-1:1); Q = Q(:,end:-1:1);
    Afun = Q * diag(f(d)) * Q';
    for k_idx = 1:length(k_list)
        k = k_list(k_idx);
        [f_idx k]
        Afuns = zeros(n,n,trials);
        for trial = 1:trials
            [U,D,jacks(k_idx,trial,f_idx)] = nystrom_jack(@(X) A*X,n,k,f);
            Afuns(:,:,trial) = U*diag(f(diag(D)))*U';
            errs(k_idx,trial,f_idx) = norm(Afun - Afuns(:,:,trial),"fro");
        end
        stds(k_idx,f_idx) = sqrt(norm(Afuns - mean(Afuns,3),"fro")^2 / trials);
    end
    
    subplot(1,2,f_idx)

    quantiles = quantile(jacks(:,:,f_idx),[0.1 0.5 0.9],2);
    plot_shaded(k_list,quantiles(:,2),quantiles(:,1),quantiles(:,3),orange,"Marker","o","MarkerSize",10)
    quantiles = quantile(errs(:,:,f_idx),[0.1 0.5 0.9],2);
    plot_shaded(k_list,quantiles(:,2),quantiles(:,1),quantiles(:,3),purple,"Marker","*","MarkerSize",10)
    plot(k_list,stds(:,f_idx),"Color",pink,"Marker","^","MarkerSize",10)
    set(gca,"YScale","log")
    xlabel("Approximation rank $k$")
    ylabel("Error quantities")
    if f_idx == 2
        title("square root","FontName","Courier New","Interpreter","TeX")
        legend({"","Jackknife","","Error","Std dev"},"Location","northeast")
    else
        title("spectral projector","FontName","Courier New","Interpreter","TeX")
    end
end

save("../data/fig_19_1.mat","stds","errs","jacks")

exportgraphics(gcf,"../figs/fig_19_1.png")
saveas(gcf,"../figs/fig_19_1.fig")