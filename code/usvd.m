function [U,S,V] = usvd(B,Bt,n,k)
% Input:  Functions B() and Bt() computing matrix products B(X) = 
%         B*X and Bt(X) = B'*X, number of rows n, and number of
%         matvecs s
% Output: Unbiased low-rank approximation Bhat to B, presented as an
%         economy-size SVD Bhat = U*S*V'

% Define test matrix 
Om = sqrt(n) * cnormc(randn(n,k)); % Matrix of random signs
Y = B(Om);              % Matvecs with B
[Q,R] = qr(Y,"econ");   % Randomized SVD
S = cnormc(inv(R'));    % Downdate randomized SVD
C = Bt(Q);              % Matvecs with B'

% Compute unbiased randomized SVD
F = (eye(k) - S*S'/k) * C' + (S .* diagprod(S,R).') * Om' / k;
[UU,S,V] = svd(F,"econ");
U = Q*UU;

end