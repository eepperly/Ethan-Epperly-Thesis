thesis_startup

load("../data/nobel.mat")
mean_age = mean(age);
centered_age = age - mean_age;

kernel = @(S,T) exp(-pdist2(S,T).^2/50);
xx = linspace(min(prize),max(prize),500)';
fhat = krr(centered_age,kernel,prize,1e-13);
plot(xx,mean_age+fhat(xx),"Color",orange); hold on
scatter(prize,age,"filled","MarkerFaceColor",blue)
axis([1900 2024 20 100])
xlabel("Year")
ylabel("Laureate Age")
set(gca,"FontSize",25)
exportgraphics(gcf,"../figs/fig_5_3_a.png")
saveas(gcf,"../figs/fig_5_3_a.fig")

figure
fhat = krr(centered_age,kernel,prize,1);
plot(xx,mean_age+fhat(xx),"Color",orange); hold on
scatter(prize,age,"filled","MarkerFaceColor",blue)
axis([1900 2024 20 100])
xlabel("Year")
set(gca,"FontSize",25)
exportgraphics(gcf,"../figs/fig_5_3_b.png")
saveas(gcf,"../figs/fig_5_3_b.fig")

figure
fhat = krr(centered_age,kernel,prize,1000);
plot(xx,mean_age+fhat(xx),"Color",orange); hold on
scatter(prize,age,"filled","MarkerFaceColor",blue)
axis([1900 2024 20 100])
xlabel("Year")
set(gca,"FontSize",25)
exportgraphics(gcf,"../figs/fig_5_3_c.png")
saveas(gcf,"../figs/fig_5_3_c.fig")