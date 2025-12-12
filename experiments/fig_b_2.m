n = 1e3;
rho = 0.5;
d = n / rho;

trials = 100;
X = zeros(n,trials);

for trial = 1:trials
    trial
    X(:,trial) = svd(random_signs(d,n)/sqrt(d));
end

%% Plot

figure('Units', 'normalized', 'Position', [0, 0.4, 0.5, 0.3]);
tiledlayout(1,3);
nexttile([1 2])

hold on

bin_edges = linspace(0, 2.5, 101);  % Define bin edges

% Loop over each column of X and plot histogram with transparency
for i = 1:size(X, 2)
    % Plot histogram of the ith column with transparency, same color, no edge
    h = histogram(X(:,i), "Normalization", "pdf", "BinEdges", bin_edges);
    h.FaceColor = darkpurple;  % Set same face color
    h.EdgeColor = 'none';  % Remove black outlines
    h.FaceAlpha = 0.02;  % Set transparency to 0.1
end

smin = 1 - sqrt(rho);
smax = 1 + sqrt(rho);
xx = linspace(smin,smax,200);
plot(xx,sqrt((smax^2-xx.^2).*(xx.^2-smin^2)) ./ (pi*rho*xx),"k--")

axis([0 2.5 0 1.5])

% Add labels and title
xlabel('$\sigma_i(\mbox{\boldmath $S$}^*\mbox{\boldmath $Q$})$')
ylabel('Probability density')

nexttile(3)
stairs(sort(1 - min(X)),(1:trials)/trials, "b", 'LineWidth', 3); hold on % Set blue solid line
stairs(sort(max(X)) - 1,(1:trials)/trials, "r--", 'LineWidth', 3) % Set red dashed line
xline(sqrt(rho), 'k:', 'LineWidth', 3) % Keeping yline black
legend({"$\eta_-$","$\eta_+$"},"Location","NorthWest","FontSize",16)
xlabel("Distortion $\eta$")
ylabel("Cumulative distribution")

% axis([0 1 0 1]) % Set axes limits

saveas(gcf,"../figs/fig_b_2.fig")
exportgraphics(gcf,"../figs/fig_b_2.png")
