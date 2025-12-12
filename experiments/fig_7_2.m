thesis_startup

f = @(S) 1 ./ (1 + sqrownorms(S).^0.5); % f(x) = 1/(1+||x||)
laplace = @(S,T) exp(-pdist2(S,T));
S = rpcholesky_rejection(laplace, @() randn(1,2), 200, 2);
fhat = kernel_interp(f,laplace,S);

[xx,yy] = meshgrid(linspace(-4,4));
z = f([xx(:) yy(:)]);
zhat = fhat([xx(:) yy(:)]);

tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile
contourf(xx,yy,reshape(z,size(xx)))
clim([0 1]); axis square

nexttile
contourf(xx,yy,reshape(zhat,size(xx)))
clim([0 1]); axis square
hold on
scatter(S(:,1),S(:,2),10,"white","filled")

cb = colorbar('Location', 'southoutside'); 
cb.Layout.Tile = 'south';

figure
colormap turbo
contourf(xx,yy,reshape(abs(z-zhat),size(xx)))
axis square
colorbar('Location', 'southoutside')

figure(1)
exportgraphics(gcf,'../figs/fig_7_2_a.png')
saveas(gcf,'../figs/fig_7_2_a.fig')
figure(2)
exportgraphics(gcf,'../figs/fig_7_2_b.png')
saveas(gcf,'../figs/fig_7_2_b.fig')