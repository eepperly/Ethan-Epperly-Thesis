n = 1e3;
m = 1e5;
trials = 100;

rhos = [0.5 n/ceil(n*log(n))];
Xs = cell(length(rhos),1);
Q = [eye(n);zeros(m-n,n)];

for j = 1:length(rhos)
    d = n / rhos(j);
    X = zeros(n,trials);
    for trial = 1:trials
        [j trial]
        X(:,trial) = svd(srtt(d,m,1)*Q);
    end
    Xs{j} = X;
end

%% Plot

figure('Units', 'normalized', 'Position', [0, 0.4, 0.5, 0.6]);
tiledlayout(2,3);

for j = 1:length(rhos)
    X = Xs{j};

    nexttile([1 2])
    hold on
    
    bin_edges = linspace(0, 2.5, 101);  % Define bin edges

    % Loop over each column of X and plot histogram with transparency
    for i = 1:size(X, 2)
        % Plot histogram of the ith column with transparency, same color, no edge
        h = histogram(X(:,i), "Normalization", "pdf", "BinEdges", bin_edges);
        h.FaceColor = blue;  % Set same face color
        h.EdgeColor = 'none';  % Remove black outlines
        h.FaceAlpha = 0.02;  % Set transparency to 0.1
    end

    smin = 1 - sqrt(rhos(j));
    smax = 1 + sqrt(rhos(j));
    xx = linspace(smin,smax,200);
    plot(xx,sqrt((smax^2-xx.^2).*(xx.^2-smin^2)) ./ (pi*rhos(j)*xx),"k--")

    axis([0 2.5 0 0.5+j])

    % Add labels and title
    if j == 2
        xlabel('$\sigma_i(\mbox{\boldmath $S$}^*\mbox{\boldmath $Q$})$') 
    end
    ylabel('Prob.\ density')

    if j == 1
        title("$d=2n$")
    else
        title("$d=\lceil n \log n\rceil$")
    end
    
    nexttile()
    stairs(sort(1 - min(X)),(1:trials)/trials, "b", 'LineWidth', 3); hold on % Set blue solid line
    stairs(sort(max(X)) - 1,(1:trials)/trials, "r--", 'LineWidth', 3) % Set red dashed line
    xline(sqrt(rhos(j)), 'k:', 'LineWidth', 3) % Keeping yline black
    if j == 2
        legend({"$\eta_-$","$\eta_+$"},"Location","SouthEast","FontSize",16)
        xlabel("Distortion $\eta$")
    end
    ylabel("Cumulative dist.")
    
    axis([0 1.5 0 1]) % Set axes limits
end

saveas(gcf,"../figs/fig_23_3.fig")
exportgraphics(gcf,"../figs/fig_23_3.png")