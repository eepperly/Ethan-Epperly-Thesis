thesis_startup

f = @(x) sin((x(:,1) + x(:,2))/10);
n = 1e4;
eye_points = 1e2;
mouth_points = 1e3;
face_points = n - 2*eye_points - mouth_points;
D = zeros(n,2);

idx = 1;
for x_shift = [-4 4]
    for i = 1:eye_points
        while true
            x = 2*rand-1;
            y = 2*rand-1;
            if x^2 + y^2 <= 1
                D(idx,1) = x + x_shift;
                D(idx,2) = y + 4;
                idx = idx + 1;
                break
            end
        end
    end
end

for x = linspace(-5,5,mouth_points)
    D(idx,1) = x;
    D(idx,2) = x^2 / 16 - 5;
    idx = idx + 1;
end

tt = linspace(0,2*pi,face_points+1);
for t = tt(1:face_points)
    D(idx,1) = 10 * cos(t);
    D(idx,2) = 10 * sin(t);
    idx = idx + 1;
end

bandwidth = 40;
kernel = @(X,Y) (1+sqrt(5)*pdist2(X,Y)/bandwidth+5*pdist2(X,Y).^2/(3*bandwidth^2)) .* exp(-sqrt(5)*pdist2(X,Y)/bandwidth);
trials = 100;
ks = 25:25:500;
rpc_errs = zeros(length(ks),trials);
unif_errs = zeros(length(ks),trials);
y = f(D);

for trial = 1:trials
    trial
    [~,Srpc] = rpcholesky(@(S) kernel(D,D(S,:)),ones(n,1),max(ks));
    Sunif = randsample(n,max(ks));
    for k_idx = 1:length(ks)
        k = k_idx;
        B = kernel(D,D(Srpc(1:k),:));
        beta = B \ y;
        rpc_errs(k_idx,trial) = norm(y - B*beta,"Inf");
        B = kernel(D,D(Sunif(1:k),:));
        beta = B \ y;
        unif_errs(k_idx,trial) = norm(y - B*beta,"Inf");
    end
end

%% Plot

figure
scatter(D(:,1), D(:,2), 40, 'k', 'filled', 'MarkerFaceAlpha', 0.1); hold on 
scatter(D(Sunif(1:100),1), D(Sunif(1:100),2), 200, hex2rgb(purple), 's','LineWidth', 2);
scatter(D(Srpc(1:100),1), D(Srpc(1:100),2), 80, hex2rgb(orange), '*','LineWidth', 1.5);
axis equal
axis off                   
set(gca, 'Color', 'none')  
set(gcf, 'Color', 'none') 
saveas(gcf,"../figs/fig_6_2_a.fig")
exportgraphics(gcf,"../figs/fig_6_2_a.png")

figure
quantiles = quantile(unif_errs,[0.1 0.5 0.9],2);
plot_shaded(ks,quantiles(:,2),quantiles(:,1),quantiles(:,3),purple,"Marker","s","MarkerFaceColor",purple,"MarkerSize",10)
quantiles = quantile(rpc_errs,[0.1 0.5 0.9],2);
plot_shaded(ks,quantiles(:,2),quantiles(:,1),quantiles(:,3),orange,"Marker","*","MarkerSize",10)
set(gca,"YScale","log")
xlabel("Number of points $k$")
ylabel("Error")
legend("","Uniform","","RPCholesky")
saveas(gcf,"../figs/fig_6_2_b.fig")
exportgraphics(gcf,"../figs/fig_6_2_b.png")