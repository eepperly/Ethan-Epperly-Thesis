n = 1e3;
rho = 0.5;
d = n / rho;

trials = 100;
zetas = [1 2 4 8 16];
Xs = cell(length(zetas),1);

for j = 1:length(zetas)
    X = zeros(n,trials);
    for trial = 1:trials
        [j trial]
        X(:,trial) = svd(full(sparse_sign(d,n,zetas(j))));
    end
    Xs{j} = X;
end

%% Plot

figure('Units', 'normalized', 'Position', [0, 0.4, 0.5, 1]);
tiledlayout(length(zetas),3);

for j = 1:length(zetas)
    X = Xs{j};

    nexttile([1 2])
    hold on
    
    bin_edges = linspace(0, 2.5, 101);  % Define bin edges

    % Loop over each column of X and plot histogram with transparency
    for i = 1:size(X, 2)
        % Plot histogram of the ith column with transparency, same color, no edge
        h = histogram(X(:,i), "Normalization", "pdf", "BinEdges", bin_edges);
        h.FaceColor = orange;  % Set same face color
        h.EdgeColor = 'none';  % Remove black outlines
        h.FaceAlpha = 0.02;  % Set transparency to 0.1
    end

    smin = 1 - sqrt(rho);
    smax = 1 + sqrt(rho);
    xx = linspace(smin,smax,200);
    plot(xx,sqrt((smax^2-xx.^2).*(xx.^2-smin^2)) ./ (pi*rho*xx),"k--")

    axis([0 2.5 0 1.5])

    % Add labels and title
    if j == length(zetas)
        xlabel('$\sigma_i(\mbox{\boldmath $SQ$})$')
    end
    ylabel('Prob.\ density')

    title(['$\' sprintf('zeta = %d',zetas(j)) '$'])
    
    nexttile(3*j)
    stairs(sort(1 - min(X)),(1:trials)/trials, "b", 'LineWidth', 3); hold on % Set blue solid line
    stairs(sort(max(X)) - 1,(1:trials)/trials, "r--", 'LineWidth', 3) % Set red dashed line
    xline(sqrt(rho), 'k:', 'LineWidth', 3) % Keeping yline black
    if j == length(zetas)
        legend({"$\eta_-$","$\eta_+$"},"Location","SouthEast","FontSize",16)
        xlabel("Distortion $\eta$")
    end
    ylabel("Cumulative dist.")
    
    axis([0.5 1.5 0 1]) % Set axes limits
end

saveas(gcf,"../figs/fig_23_6.fig")
exportgraphics(gcf,"../figs/fig_23_6.png")