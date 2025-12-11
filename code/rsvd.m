function [U,S,V] = rsvd(B,Bt,n,k)
% Input:  Functions B() and Bt() computing matrix products B(X) = 
%         B*X and Bt(X) = B'*X, number of columns n, and rank k
% Output: Low-rank approximation Bhat to B, presented as an
%         economy-size SVD Bhat = U*S*V'

Om = randn(n,k);            % Gaussian random test matrix
Y = B(Om);                  % Matvecs with B
[Q,~] = qr(Y,"econ");       % Orthogonalize
C = Bt(Q);                  % Matvecs with B'
[UU,S,V] = svd(C',"econ");  % SVD of factor matrix
U = Q*UU;                   % Get left singular vectors

end