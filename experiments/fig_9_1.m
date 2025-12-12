thesis_startup

n = 50;
x = linspace(0,2,2*n);
x1 = x(1:n)'; x2 = x(n+1:end)';
y = linspace(0,1,n)';
X = [kron(x1,ones(n,1)) kron(ones(n,1),y)];
Y = [kron(x2,ones(n,1)) kron(ones(n,1),y)];
B = 1 ./ pdist2(X,Y);

k_list = 5:5:50;
trials = 100;

[QQ,FF] = greedy_qr(B,max(k_list));

rpqr_errs = zeros(length(k_list),trials);
cpqr_errs = zeros(length(k_list),1);
sp_errs = zeros(length(k_list),trials);

for k_idx = 1:length(k_list)
    k = k_list(k_idx)
    cpqr_errs(k_idx) = norm(B - QQ(:,1:k)*FF(:,1:k)',"fro") / norm(B,"fro");
    for trial = 1:trials
        [Q,F] = acc_rpqr(B,k,ceil(k/2));
        rpqr_errs(k_idx,trial) = norm(B - Q*F',"fro")/norm(B,"fro");

        St = sparse_sign(2*k,size(B,1),4);
        StB = St*B;
        [~,~,s] = lu(StB(1:k,:)',"vector");
        [Q,~] = qr(B(:,s(1:k)),"econ");
        sp_errs(k_idx,trial) = norm(B - Q*(Q'*B),"fro")/norm(B,"fro");
    end
end 

%% Plot

figure(1)

quantiles = quantile(sp_errs,[0.1 0.5 0.9],2);
plot_shaded(k_list,quantiles(:,2),quantiles(:,1),quantiles(:,3),purple,"Marker","*","MarkerSize",10)
plot(k_list,cpqr_errs,"Color",blue,"Marker","x","MarkerSize",10)
quantiles = quantile(rpqr_errs,[0.1 0.5 0.9],2);
plot_shaded(k_list,quantiles(:,2),quantiles(:,1),quantiles(:,3),orange,"Marker","o","MarkerSize",10)

set(gca,"YScale","log")

xlabel("Approximation rank $k$")
ylabel("Relative Frobenius-norm error")

legend({"","Sketchy pivoting","Greedy pivoting (CPQR)", "", "RPQR"},"Location","SouthWest")

exportgraphics(gcf,"../figs/fig_9_1.png")
saveas(gcf,"../figs/fig_9_1.fig")