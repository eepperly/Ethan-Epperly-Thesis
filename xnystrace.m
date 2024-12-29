function [tr,est] = xnystrace(A,n,s)
% Input:  Function A() computing matrix products A(X) = A*X, number
%         of rows n, and number of matvecs s
% Output: Estimate tr of trace(A), estimate est of the error
%         abs(tr - trace(A))

% Single-pass Nystrom approximation
Om = random_signs(n,s);         % Matrix of random signs
Y = A(Om);                      % Matrix product Y = A*Om
mu = eps*norm(Y,"fro")/sqrt(n); % Compute shift
Y = Y + mu * Om;                % Apply shift to Y
H = Om'*Y;
R = chol((H+H')/2);             % Explicitly symmetrize H to be safe
F = Y/R;                        % Triangular substitution

% Apply Nystrom downdating
Z = (F/R') .* (sqrownorms(inv(R)) .^ (-1/2))';

% Compute estimator and estimate
tr_vec = norm(F,"fro")^2 * ones(s,1) - sqcolnorms(Z) ...
    + abs(diagprod(Om,Z)) .^ 2 - mu * n * ones(s,1);
tr = mean(tr_vec);           % Trace estimate
est = std(tr_vec) / sqrt(s); % Error estimate

end

