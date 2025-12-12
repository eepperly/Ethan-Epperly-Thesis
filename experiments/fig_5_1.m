thesis_startup

kernel = @(S,T) cosh(min(S,T')) .* cosh(1-max(S,T')) / sinh(1);
x = linspace(0,1,500);
subplot(1,2,1);
plot(x,kernel(x',0.25),"Color",orange); hold on
plot(x,kernel(x',0.5),"Color",blue)
subplot(1,2,2)
[xx,yy] = meshgrid(x); zz = kernel(x',x');
contourf(xx,yy,zz)
colorbar

fig = gcf;  % Get the current figure handle
position = fig.Position;  % Get the current position [left, bottom, width, height]
position(3) = 2 * position(3);  % Double the width (third element)
fig.Position = position;  % Set the new position
exportgraphics(fig,"../figs/fig_5_1.png")
saveas(fig,"../figs/fig_5_1.fig")