thesis_startup

n = 50;
x = linspace(0,2,2*n);
x1 = x(1:n)'; x2 = x(n+1:end)';
y = linspace(0,1,n)';
X = [kron(x1,ones(n,1)) kron(ones(n,1),y)];
Y = [kron(x2,ones(n,1)) kron(ones(n,1),y)];
B = 1 ./ pdist2(X,Y);

ks = 100:100:2500;

rpqr_times = zeros(length(ks),1);
acc_rpqr_times = zeros(length(ks),1);

for k_idx = 1:length(ks)
    k = ks(k_idx)
    tic; rpqr(B,k); rpqr_times(k_idx) = toc;
    rpqr_times(k_idx)
    tic; acc_rpqr(B,k,min(k/2,200)); acc_rpqr_times(k_idx) = toc;
    acc_rpqr_times(k_idx)
end

tic; [Q,R,P] = qr(B); cpqr_time = toc;

%% Plot

figure

plot(ks,rpqr_times,"Color",blue,"Marker","x"); hold on
plot(ks,acc_rpqr_times,"Color",orange,"Marker","o");
scatter(n^2,cpqr_time,100,hex2rgb(purple),"Marker","*","LineWidth",2)
set(gca,"YScale","log")
xlabel("Approximation rank $k$")
ylabel("Time (sec)")
legend({"RPQR","Acc.\ RPQR","Greedy (CPQR)"},"Location","SouthEast")

exportgraphics(gcf,"../figs/fig_9_2.png")
saveas(gcf,"../figs/fig_9_2.fig")