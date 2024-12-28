function [tr, est] = xtrace(B,n,s)
% Input:  Function B() computing matrix products B(X) = B*X, number
%         of rows n, and number of matvecs s
% Output: Estimate tr of trace(B), estimate est of the error
%         abs(tr - trace(B))

% Define test matrix 
k = floor(s/2);         % Approximation rank is s/2
Om = random_signs(n,k); % Matrix of random signs

% Randomized SVD and downdate
Y = B(Om);            % Collect matvecs
[Q,R] = qr(Y,"econ"); % Randomized SVD
S = cnormc(inv(R'));  % cnormc normalizes the columns

% Compute other necessary matrices
Z = B(Q);  % Collect matvecs
H = Q'*Z;
W = Q'*Om;
T = Z'*Om;
X = W - diagprod(S,W) .* S;

% Compute estimator, output
tr_vec = trace(H) * ones(k,1) - diagprod(H*S,S) - diagprod(X,T)...
            + diagprod(H*X,X) + diagprod(S, W) .* diagprod(R, S);
tr = mean(tr_vec);           % Trace estimate
est = std(tr_vec) / sqrt(k); % Error estimate

end
