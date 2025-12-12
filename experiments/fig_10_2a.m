thesis_startup

n = 1000;
x = exp(pi*1i*(0:(n-1))/n);
y = exp(pi*1i*(n:(2*n-1))/n);
B = 1 ./ (x-y.');
[Z,D,V] = svd(B);
Bfro = norm(B,"fro");

ks = 5:5:80;
num_trials = 100;

rpcur2_errs = zeros(length(ks),num_trials);
rpcurlev_errs = zeros(length(ks),num_trials);
rpcurlev_noover_errs = zeros(length(ks),num_trials);
mdcur_errs = zeros(length(ks),num_trials);
rpcur2_num_errs = zeros(length(ks),num_trials);
rpcurlev_num_errs = zeros(length(ks),num_trials);
rpcurlev_noover_num_errs = zeros(length(ks),num_trials);
mdcur_num_errs = zeros(length(ks),num_trials);
best = zeros(length(ks));

warning('off', 'all')
for k_idx = 1:length(ks)
    k = ks(k_idx)
    l = min(ceil(1.5*log(k)*k),n);
    kk = ceil(k/(2*log(k)));
    for trial = 1:num_trials
        [S,T,U,G,P] = rpcur2(B,k,k);
        rpcur2_errs(k_idx,trial) =...
            norm(B - B(:,S)*(P\(G*B(T,:))),"fro") / norm(B,"fro");
        rpcur2_num_errs(k_idx,trial) =...
            norm(B - B(:,S)*U*B(T,:),"fro") / norm(B,"fro");
        [S,T,U,G,P] = rpcur_lev(B,k,l);
        rpcurlev_errs(k_idx,trial) =...
            norm(B - B(:,S)*(P\(G*B(T,:))),"fro") / Bfro;
        rpcurlev_num_errs(k_idx,trial) =...
            norm(B - B(:,S)*U*B(T,:),"fro") / norm(B,"fro");
        [S,T,U,G,P] = rpcur_lev(B,k,k);
        rpcurlev_noover_errs(k_idx,trial) =...
            norm(B - B(:,S)*(P\(G*B(T,:))),"fro") / Bfro;
        rpcurlev_noover_num_errs(k_idx,trial) =...
            norm(B - B(:,S)*U*B(T,:),"fro") / norm(B,"fro");
        [S,T,U,G,P] = md_cur(B,k,l,Z(:,1:kk),V(:,1:kk));
        mdcur_errs(k_idx,trial) =...
            norm(B - B(:,S)*(P\(G*B(T,:))),"fro") / Bfro;
        mdcur_num_errs(k_idx,trial) =...
            norm(B - B(:,S)*U*B(T,:),"fro") / norm(B,"fro");
    end
    best(k_idx) = norm(B - Z(:,1:k)*D(1:k,1:k)*V(:,1:k)',"fro") / Bfro;
end
warning('on', 'all')

quantiles = quantile(rpcurlev_noover_errs,[0.1 0.5 0.9],2);
plot_shaded(ks,quantiles(:,2),quantiles(:,1),quantiles(:,3),yellow,"Marker","s","MarkerSize",10,"MarkerFaceColor",yellow)
quantiles = quantile(rpcurlev_noover_num_errs,[0.1 0.5 0.9],2);
plot_shaded(ks,quantiles(:,2),quantiles(:,1),quantiles(:,3),darkyellow,"Marker","s","MarkerSize",10,"LineStyle","--","MarkerFaceColor",darkyellow)

quantiles = quantile(mdcur_errs,[0.1 0.5 0.9],2);
plot_shaded(ks,quantiles(:,2),quantiles(:,1),quantiles(:,3),pink,"Marker","x","MarkerSize",10)
quantiles = quantile(mdcur_num_errs,[0.1 0.5 0.9],2);
plot_shaded(ks,quantiles(:,2),quantiles(:,1),quantiles(:,3),darkpink,"Marker","x","MarkerSize",10,"LineStyle","--")

quantiles = quantile(rpcurlev_errs,[0.1 0.5 0.9],2);
plot_shaded(ks,quantiles(:,2),quantiles(:,1),quantiles(:,3),blue,"Marker","o","MarkerSize",10,"MarkerFaceColor",blue)
quantiles = quantile(rpcurlev_num_errs,[0.1 0.5 0.9],2);
plot_shaded(ks,quantiles(:,2),quantiles(:,1),quantiles(:,3),darkblue,"Marker","o","MarkerSize",10,"LineStyle","--","MarkerFaceColor",darkblue)

quantiles = quantile(rpcur2_errs,[0.1 0.5 0.9],2);
plot_shaded(ks,quantiles(:,2),quantiles(:,1),quantiles(:,3),orange,"Marker","*","MarkerSize",10)
quantiles = quantile(rpcur2_num_errs,[0.1 0.5 0.9],2);
plot_shaded(ks,quantiles(:,2),quantiles(:,1),quantiles(:,3),darkorange,"Marker","*","MarkerSize",10,"LineStyle","--")

plot(ks,best,"k:")
set(gca,"YScale","log")

axis([-Inf Inf 1e-15 1e0])

xlabel("Rank $k$")
ylabel("Relative error $\|\mbox{\boldmath $B$} - \mbox{\boldmath $\widehat{B}$}\|_{\rm F}/\|\mbox{\boldmath $B$}\|_{\rm F}$")

legend({"","RPCURLev ($\ell=k$)","","",...
    "","MDCUR","","",...
    "","RPCURLev ($\ell>k$)","","",...
    "","RPCUR2","","",...
    "Optimal"}, "Location", "southwest","FontSize",16)

exportgraphics(gcf,"../figs/fig_10_2_a.png")
saveas(gcf,"../figs/fig_10_2_a.fig")
