thesis_startup

load("../data/SUSY.mat")

shuffled = randperm(size(A,1));
n = 1e5;
D = A(shuffled(1:n),:);
y = b(shuffled(1:n),:);
Dtest = A(shuffled(n+1:(2*n)),:);
ytest = b(shuffled(n+1:(2*n)));

% Standardize
meanD = mean(D); stdD = std(D);
D = (D - meanD) ./ stdD; Dtest = (Dtest - meanD) ./ stdD;

kernel = @(X,Y) exp(-pdist2(X,Y)/8);

ks = 25:25:1000;
active_rpcs_accs = zeros(size(ks));
active_unif_accs = zeros(size(ks));
restricted_rpcs_accs = zeros(size(ks));
restricted_unif_accs = zeros(size(ks));
[~,~,Srpc] = rpcholesky_active_krr(y,kernel,D,0,max(ks));
Sunif = randsample(n,max(ks));

for k_idx = 1:length(ks)
    k = ks(k_idx)
    beta = kernel(D(Srpc(1:k),:),D(Srpc(1:k),:)) \ y(Srpc(1:k));
    active_rpcs_accs(k_idx) = 1 - mean(ytest == (kernel(Dtest,D(Srpc(1:k),:)) * beta > 0.5));
    beta = kernel(D,D(Srpc(1:k),:)) \ y;
    restricted_rpcs_accs(k_idx) = 1 - mean(ytest == (kernel(Dtest,D(Srpc(1:k),:)) * beta > 0.5));
    beta = kernel(D(Sunif(1:k),:),D(Sunif(1:k),:)) \ y(Sunif(1:k));
    active_unif_accs(k_idx) = 1 - mean(ytest == (kernel(Dtest,D(Sunif(1:k),:)) * beta > 0.5));
    beta = kernel(D,D(Sunif(1:k),:)) \ y;
    restricted_unif_accs(k_idx) = 1 - mean(ytest == (kernel(Dtest,D(Sunif(1:k),:)) * beta > 0.5));
end

beta = kernel(D,D(Srpc(1:k),:)) \ y;
restricted_acc_rpc = 1 - mean(ytest == (kernel(Dtest,D(Srpc(1:k),:)) * beta > 0.5));

%% Plot

figure
plot(ks, active_unif_accs, "Color", purple, "Marker", "s", "MarkerFaceColor", purple, "MarkerSize", 10); hold on
plot(ks, active_rpcs_accs, "Color", orange, "Marker", "*", "MarkerSize", 10)
plot(ks, restricted_unif_accs, "--", "Color", darkpurple, "Marker", "s", "MarkerFaceColor", darkpurple, "MarkerSize", 10)
plot(ks, restricted_rpcs_accs, "--", "Color", darkorange, "Marker", "*", "MarkerSize", 10)

axis([0 1000 0.19 0.38])
xlabel("Subset size $k$")
ylabel("Test misclassification rate")
legend({"Uniform (active)", "RPCholesky (active)", "Uniform (restricted)", "RPCholesky (restricted)"},"Location","northeast")

saveas(gcf, "../figs/fig_6_3.fig")
exportgraphics(gcf, "../figs/fig_6_3.png")