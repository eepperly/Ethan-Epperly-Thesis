function d = udiagpp(B,Bt,n,s)
% Input:  Functions B() and Bt() computing matrix products B(X) = 
%         B*X and Bt(X) = B'*X, number of rows n, and number of
%         matvecs s
% Output: Estimate d of the diagonal of B

% Randomized SVD
k = floor(s/3);         % Approximation rank is s/3
Om = random_signs(n,k); % Matrix of random signs
BOm = B(Om);            % Matvecs with B
[Q,~] = qr(BOm,"econ"); % Randomized SVD

% Unbiased Diag++ estimator
BtQ = Bt(Q);            % Matvecs with B'
k = s - 2*k;            % Remaining matvecs
Ga = random_signs(n,k); % More random signs
BGa = B(Ga);            % Matvecs with B
d = diagprod(Q',BtQ') + mean((BGa - Q*(Q'*BGa)) .* Ga,2);

end