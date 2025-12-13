thesis_startup

iter = 50;

[B,c,x,r,s,V] = random_ls_problem(4000,50,1e12,1e-4);
St = sparse_sign(1000,4000,4);

colors = {yellow,purple,blue,pink,orange};
styles = {"-","--",":","-.","-"};
markers = {"s","*","x","^","o"};

for i = 1:5
    if i == 1
        X = sketch_precondition_zero(B,c,St,iter);
    elseif i == 2
        X = sketch_descend(B,c,St,1,0,iter);
    elseif i == 3
        beta = size(B,2) / size(St,1);
        alpha = (1-beta)^2;
        X = sketch_descend(B,c,St,alpha,beta,iter);
    elseif i == 4
        X = sketch_precondition(B,c,St,iter);
    else
        X = spir(B,c,St,iter/2,iter/2);
    end

    figure(1)
    title("forward error","FontName","Courier New","Interpreter","TeX")
    semilogy(0:iter,sqrt(sqcolnorms(X - x*ones(1,iter+1)))/norm(x),...
        styles{i},"Color",colors{i},"LineWidth",3, 'MarkerIndices', 1:5:(iter+1),"Marker",markers{i},"MarkerSize",15); hold on
    if i == 5
        yline(norm(x - B\c)/norm(x),"k:","LineWidth",4)
    end
    ylabel("Relative forward error")
    xlabel("Iteration")
    axis([0 50 1e0 1e10])
    
    figure(2)
    title("residual error","FontName","Courier New","Interpreter","TeX")
    semilogy(0:iter,sqrt(sqcolnorms(B*(X - x*ones(1,iter+1))))/norm(B*x),...
        styles{i},"Color",colors{i},"LineWidth",3, 'MarkerIndices', 1:5:(iter+1),"Marker",markers{i},"MarkerSize",15); hold on
    if i == 5
        yline(norm(B*(x - B\c))/norm(B*x),"k:","LineWidth",4)
    end
    ylabel("Relative residual error")
    xlabel("Iteration")
    axis([0 50 1e-10 1e0])

    if i == 5
        legend({"S\&P ($\mbox{\boldmath $x$}_0 = \mbox{\boldmath $0$}$)", "S\&D (trivial)", "S\&D (opt)", "S\&P", "SPIR", "Direct"},"Location","northeast")
    end
    
    figure(3)
    title("backward error","FontName","Courier New","Interpreter","TeX")
    bes = zeros(iter+1,1);
    for j = 1:(iter+1)
        bes(j) = backerr_est(B,c,X(:,j),s,V) / norm(B,"fro");
    end
    semilogy(0:iter,bes,styles{i},"Color",colors{i},"LineWidth",3, 'MarkerIndices', 1:5:(iter+1),"Marker",markers{i},"MarkerSize",15); hold on
    if i == 5
        yline(backerr_est(B,c,B\c,s,V)/norm(B,"fro"),"k:","LineWidth",4)
    end
    ylabel("Backward error estimate")
    xlabel("Iteration")
    axis([0 50 1e-17 1e-10])
end

figure(1)
saveas(gcf,"../figs/fig_22_1_a.fig")
exportgraphics(gcf,"../figs/fig_22_1_a.png")

figure(2)
saveas(gcf,"../figs/fig_22_1_b.fig")
exportgraphics(gcf,"../figs/fig_22_1_b.png")

figure(3)
saveas(gcf,"../figs/fig_22_1_c.fig")
exportgraphics(gcf,"../figs/fig_22_1_c.png")