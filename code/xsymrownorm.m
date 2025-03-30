function srn = xsymrownorm(A,n,s)
% Input:  Function A() computing matrix products A(X) = A*X, number
%         of rows n, and number of matvecs s
% Output: Estimate srn of the squared row norms of A

% Randomized SVD and downdating
k = floor(s/2);         % Approximation rank is s/2
Om = random_signs(n,k); % Matrix of random signs
Y = A(Om);              % Matvecs 
[Q,R] = qr(Y,"econ");   % Orthogonalize
S = cnormc(inv(R'));    % cnormc normalizes columns

% Compute other necessary matrices
Z = A(Q);               % Matvecs
W = Q'*Om;
X = W - S .* diagprod(S,W).';

% Form estimate
srn = sqrownorms(Z) + (-sqrownorms(Z*S) + sqrownorms(Y - Z*X))/k;

end