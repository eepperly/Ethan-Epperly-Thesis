function [U,S,V,est] = rsvd_errest(B,Bt,n,k)
% Input:  Functions B() and Bt() computing matrix products B(X) = 
%         B*X and Bt(X) = B'*X, number of columns n, and rank k
% Output: Low-rank approximation Bhat to B, presented as an
%         economy-size SVD Bhat = U*S*V', error estimate

Om = random_signs(n,k);                  % Matrix of random signs
Y = B(Om);                               % Matvecs with B
[Q,R] = qr(Y,"econ");                    % Orthogonalize
C = Bt(Q);                               % Matvecs with B'
[UU,S,V] = svd(C',"econ");               % SVD of factor matrix
U = Q*UU;                                % Get left singular vectors
est = sqrt(mean(1./sqrownorms(inv(R)))); % Error estimate

end