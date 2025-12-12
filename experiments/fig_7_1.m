thesis_startup

laplace = @(S,T) exp(-pdist2(S,T)/10);
S = rpcholesky_rejection(laplace,@() rand(1,2),50,2);
scatter(S(:,1),S(:,2),100,"filled","MarkerFaceColor",orange)
axis off; figure
scatter(rand(1,50),rand(1,50),100,"filled","MarkerFaceColor",blue)
axis off

saveas(1,'../figs/fig_7_1_a.png')
saveas(1,'../figs/fig_7_1_a.fig')
saveas(2,'../figs/fig_7_1_b.png')
saveas(2,'../figs/fig_7_1_b.fig')