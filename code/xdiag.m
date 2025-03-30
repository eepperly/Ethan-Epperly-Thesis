function d = xdiag(B,Bt,n,s)
% Input:  Functions B() and Bt() computing matrix products B(X) = 
%         B*X and Bt(X) = B'*X, number of rows n, and number of
%         matvecs s
% Output: Estimate d of the diagonal of B

% Randomized SVD
k = floor(s/2);         % Approximation rank is s/2
Om = random_signs(n,k); % Matrix of random signs
Y = B(Om);              % Matvecs with B
[Q,R] = qr(Y,"econ");   % Randomized SVD
S = cnormc(inv(R'));    % Downdate RSVD, cnormc normalizes cols

% Form estimator
Z = Bt(Q);              % Matvecs with B'
d = diagprod(Q',Z') + mean((Q*S)...                   % XDiag
    .* (-conj(Z*S) + conj(Om) .* diagprod(S,R).'),2);

end