thesis_startup

n = 1000;
x = exp(pi*1i*(0:(n-1))/n);
y = exp(pi*1i*(n:(2*n-1))/n);
B = 1 ./ (x-y.');

k = 30; l = ceil(1.5*k*log(k));

figure('Units', 'normalized', 'Position', [0.1 0.1 0.4 0.8]);
t = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:2
    nexttile
    if i == 1
        [S,T] = rpcur2(B,k,k);
        color = orange;
    else
        [U,~,V] = svd(B);
        kk = ceil(k/(2*log(k)));
        [S,T] = md_cur(B,k,k,U(:,1:kk),V(:,1:kk));
        color = pink;
    end

    % scatter(real(x),imag(x),20,"filled","MarkerFaceColor",black,"MarkerFaceAlpha",0.05); hold on
    scatter(real(y),imag(y),20,"filled","MarkerFaceColor",black,"MarkerFaceAlpha",0.05); hold on
    % scatter(real(x(T)),imag(x(T)),400,"filled","Color",color,"MarkerFaceColor",color,"MarkerFaceAlpha",0.5,"Marker","o")
    scatter(real(y(S)),imag(y(S)),400,"filled","Color",color,"MarkerFaceColor",color,"MarkerFaceAlpha",0.5,"Marker","s","MarkerEdgeColor",black,"LineWidth",1.5)

    axis equal
    axis([-1.1 1.1 -1.1 0])
    axis off

    if i == 1
        title("Random pivoting","FontSize",30)
    else
        title("Mahoney--Drineas","FontSize",30)
    end
end

exportgraphics(gcf,"../figs/fig_10_2_b.png")
saveas(gcf,"../figs/fig_10_2_b.fig")
