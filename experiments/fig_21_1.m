thesis_startup

iter = 50;

[B,c,x,r,s,V] = random_ls_problem(4000,50,1e12,1e-4);
St = sparsesign(100,4000,4);
X = sketch_precondition_zero(B,c,St,iter);

figure("Position", [100, 100, 1350, 350])

subplot(1,3,1)
semilogy(0:iter,sqrt(sqcolnorms(X - x*ones(1,iter+1)))/norm(x),...
    "-.","Color",yellow,"LineWidth",3)
yline(norm(x - B\c)/norm(x),"k:","LineWidth",4)
ylabel("Relative forward error")
xlabel("Iteration")
axis([0 50 1e0 1e10])

subplot(1,3,2)
semilogy(0:iter,sqrt(sqcolnorms(B*(X - x*ones(1,iter+1))))/norm(B*x),...
    "-.","Color",yellow,"LineWidth",3)
yline(norm(B*(x - B\c))/norm(B*x),"k:","LineWidth",4)
ylabel("Relative residual error")
xlabel("Iteration")
axis([0 50 1e-10 1e0])

subplot(1,3,3)
bes = zeros(iter+1,1);
for i = 1:(iter+1)
    bes(i) = backerr_est(B,c,X(:,i),s,V) / norm(B,"fro");
end
semilogy(0:iter,bes,"-.","Color",yellow,"LineWidth",3)
yline(backerr_est(B,c,B\c,s,V)/norm(B,"fro"),"k:","LineWidth",4)
ylabel("Backward error estimate")
xlabel("Iteration")
axis([0 50 1e-17 1e-10])
legend({"S\&P ($\mbox{\boldmath $x$}_0 = \mbox{\boldmath $0$}$)", "Direct"},"Location","east")

saveas(gcf,"../figs/fig_21_1.fig")
exportgraphics(gcf,"../figs/fig_21_1.png")