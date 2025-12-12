thesis_startup

kernel = @(S,T) exp(-pdist2(S,T).^2/2);
[S,R] = rpcholesky_rejection(kernel,@() randn(1,2),30,2);

figure
nplt = 100;
xx = linspace(-4,4,nplt);
yy = linspace(-4,4,nplt);
[xx,yy] = meshgrid(xx,yy);
zz = zeros(size(xx));
for i = 1:nplt
    for j = 1:nplt
        x = xx(i,j);
        y = yy(i,j);
        zz(i,j) = kernel([x y],[x y]) - norm(R'\kernel(S,[x y]))^2;
    end
end
contourf(xx,yy,zz,20); hold on
scatter(S(:,1),S(:,2),100,'filled','w')
% colormap('cool')
axis square
clim([0 1])
colorbar

saveas(gcf,"../figs/fig_7_3_a.png")
saveas(gcf,"../figs/fig_7_3_a.fig")