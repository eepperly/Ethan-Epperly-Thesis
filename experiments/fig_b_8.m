thesis_startup
n = 100;

% Get matrix metadata
index = ssget;

% Filter for 10,000 <= nrows <= 100,000
row_limits = [1e4, 1e5];
selected = find([index.nrows] >= row_limits(1) & [index.nrows] <= row_limits(2));

% Limit to at most 1000 matrices
max_matrices = 1000;
selected = selected(1:min(end, max_matrices));

fprintf('Processing up to %d matrices with between %d and %d rows.\n', ...
    length(selected), row_limits(1), row_limits(2));

distortions = zeros(length(selected),1);
s_vals = zeros(n,length(selected));

% Loop through selected matrices
for k = 1:length(selected)
    idx = selected(k);

    % Download and load the matrix
    Prob = ssget(idx);
    A = Prob.A;

    % Compute sketch and SVD
    U = rsvd_errest(@(X) A*X, @(X) A'*X, size(A,2), n);
    St = sparse_sign(2*n, size(U,1), 4);
    s_vals(:,k) = svd(St * U);

    distortions(k) = max([max(s_vals(:,k)) - 1, 1 - min(s_vals(:,k))]);
    [k distortions(k)]
end

% Plot histogram of distortions
figure;
histogram(distortions, 'FaceColor', orange, 'EdgeColor', 'none');
hold on;
xline(1/sqrt(2), 'k--', 'LineWidth', 1.5);  % black dashed line
xlabel('Distortion $\eta$', 'Interpreter', 'latex');
ylabel('Count');
grid on;