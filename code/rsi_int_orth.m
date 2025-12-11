function [U,S,V] = rsi_int_orth(B,Bt,m,n,k,q)
% Input:  Functions B() and Bt() computing matrix products B(X) = 
%         B*X and Bt(X) = B'*X, dimensions m and n, rank k, and 
%         number of subspace iteration steps q
% Output: Low-rank approximation Bhat to B, presented as an
%         economy-size SVD Bhat = U*S*V'

if mod(q,2) == 1            % Odd case
    Om = randn(m,k);        % Gaussian random test matrix with m rows
    Om = Bt(Om);            % First subspace iteration step
else                        % Even case
    Om = randn(n,k);        % Gaussian random test matrix with n rows
end

for i = 1:floor(q/2)
    [Om,~] = qr(Bt(B(Om)),"econ"); % Intermediate orthogonalization
end

Y = B(Om);                  % Matvecs with B
[Q,~] = qr(Y,"econ");       % Orthogonalize
C = Bt(Q);                  % Matvecs with B'
[UU,S,V] = svd(C',"econ");  % SVD of factor matrix
U = Q*UU;                   % Get left singular vectors

end