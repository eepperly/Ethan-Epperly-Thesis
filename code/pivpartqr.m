function [Q,F] = pivpartqr(B,S)
% Input:  Matrix B and set of pivots S
% Output: Factors Q and F defining a rank-k approximation Bhat = Q*F'

[m,n] = size(B);
k = length(S);                      % Rank = number of pivots
Q = zeros(m,k);                     % Left factor (orthonormal cols)
F = zeros(n,k);                     % Right factor

for i = 1:k
    s = S(i);                       % Get selected pivot
    Q(:,i) = B(:,s) / norm(B(:,s)); % Normalize column
    F(:,i) = B' * Q(:,i);           % Column of factor matrix
    B = B - Q(:,i) * F(:,i)';       % Modified Gram-Schmidt
end

end