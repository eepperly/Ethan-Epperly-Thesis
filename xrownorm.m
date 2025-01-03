function srn = xrownorm(B,Bt,n,s)
% Input:  Functions B() and Bt() computing matrix products B(X) = 
%         B*X and Bt(X) = B'*X, number of rows n, and number of
%         matvecs s
% Output: Estimate srn of the squared row norms of B

% Randomized SVD and downdating
k = floor(s/3);         % Approximation rank is s/3
Om = random_signs(n,k); % Matrix of random signs
G = B(Om);              % Matvecs with B
Y = Bt(G);              % Matvecs with B'
[Q,R] = qr(Y,"econ");   % Orthogonalize
S = cnormc(inv(R'));    % cnormc normalizes columns

% Compute other necessary matrices
Z = B(Q);                   % Matvecs with B
W = Q'*Om;
X = W - diagprod(W,S) .* S;

% Form estimate
srn = sqrownorms(Z) + (-sqrownorms(Z*S) + sqrownorms(G - Z*X))/k;

end