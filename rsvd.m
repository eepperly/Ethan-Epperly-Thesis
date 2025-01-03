function [U,S,V] = rsvd(B,Bt,n,s)
% Input:  Functions B() and Bt() computing matrix products B(X) = 
%         B*X and Bt(X) = B'*X, number of rows n, and number of
%         matvecs s
% Output: Low-rank approximation Bhat to B, presented as an
%         economy-size SVD Bhat = U*S*V'

Om = randn(n,s);            % Normal random test matrix
Y = B(Om);                  % Matvecs with B
Q = orth(Y);                % Orthogonalize
C = Bt(Q);                  % Matvecs with B'
[UU,S,V] = svd(C',"econ");  % SVD of factor matrix
U = Q*UU;                   % Get left singular vectors

end