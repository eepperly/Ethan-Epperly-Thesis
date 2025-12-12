thesis_startup
n = 1000; s = 100;
ss = 5:5:100;
trials = 100;

powers = linspace(0,3,7);
base_color = [100, 143, 255]/255;
dark_color = base_color * 0.5;
light_color = base_color + (1 - base_color) * 0.5; 
custom_colormap1 = [linspace(dark_color(1), light_color(1), length(powers))', ...
                   linspace(dark_color(2), light_color(2), length(powers))', ...
                   linspace(dark_color(3), light_color(3), length(powers))'];

base_color = [220,38,127]/255;
dark_color = base_color * 0.5;
light_color = base_color + (1 - base_color) * 0.5; 
custom_colormap2 = [linspace(dark_color(1), light_color(1), length(powers))', ...
                   linspace(dark_color(2), light_color(2), length(powers))', ...
                   linspace(dark_color(3), light_color(3), length(powers))'];


figure("Position", [100, 100, 900, 400]);
for idx = 1:length(powers)
    power = powers(idx)
    A = rand_with_evals((1:n) .^ (-power));
    hutchpp_errs = zeros(size(ss));
    xtrace_errs = zeros(size(ss));
    for i = 1:length(ss)
        s = ss(i);
        for trial = 1:trials
            hutchpp_errs(i) = hutchpp_errs(i) + abs(trace(A) - hutchpp(@(X) A*X,n,s))^2/trace(A)^2/trials;
            xtrace_errs(i) = xtrace_errs(i) + abs(trace(A) - xtrace(@(X) A*X,n,s))^2/trace(A)^2/trials;
        end
    end
    subplot(1,2,1)
    loglog(ss,sqrt(hutchpp_errs),"Color",custom_colormap1(idx,:)); hold on
    subplot(1,2,2)
    loglog(ss,sqrt(xtrace_errs),"Color",custom_colormap2(idx,:)); hold on
end

for i = 1:2
    subplot(1,2,i)
    axis([-Inf Inf 5e-6 1e0])
    plot(ss, 1./ss, "k--")
    xlabel("Number of matvecs $s$")
    if i == 1
        ylabel("Root-mean-square error")
    end
end

exportgraphics(gcf,'../figs/fig_15_1.png')
saveas(gcf,'../figs/fig_15_1.fig')