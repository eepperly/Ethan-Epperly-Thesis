thesis_startup

m_max = 1000;
s = 1;
xx = linspace(0,1,200);
f = zeros(1,200); 
g = zeros(1,200);
for m = 1:m_max
    f = f + cos(2*pi*m*xx)/m^(s+0.6) + sin(2*pi*m*xx)/m^(s+0.6);
    g = g + randn*cos(2*pi*m*xx)/m^s + randn*sin(2*pi*m*xx)/m^s;
end
plot(xx,f,"Color",orange)
figure
plot(xx,g,"Color",blue)

figure(1)
exportgraphics(gcf,"../figs/fig_5_2_a.png")
saveas(gcf,"../figs/fig_5_2_a.fig")
figure(2)
exportgraphics(gcf,"../figs/fig_5_2_b.png")
saveas(gcf,"../figs/fig_5_2_b.fig")