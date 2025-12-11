function [U,D,V,est] = rsi_errest(B,Bt,n,k,q)
% Input:  Functions B() and Bt() computing matrix products B(X) = 
%         B*X and Bt(X) = B'*X, number of columns n, rank k, and 
%         number of subspace iteration steps q
% Output: Low-rank approximation Bhat to B, presented as an
%         economy-size SVD Bhat = U*S*V' and error estimate est

assert(mod(q,2) == 0);      % Only implemented for even steps
Om = randn(n,k);            % Gaussian random test matrix with n rows

Z = B(Om); Y = Z;           % First subspace iteration step
for i = 1:(q/2-1)
    Y = B(Bt(Y));           % Two subspace iteration steps
end

[Q,R] = qr(Y,"econ");       % Orthogonalize
C = Bt(Q);                  % Matvecs with B'
[UU,D,V] = svd(C',"econ");  % SVD of factor matrix
U = Q*UU;                   % Get left singular vectors

S = cnormc(inv(R'));        % Downdating matrix
QtZ = Q'*Z;
est = norm(Z-Q*QtZ+Q*S.*diagprod(S,QtZ).',"fro")/sqrt(k);

end