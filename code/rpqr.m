function [Q,F,S] = rpqr(B,k)
% Input:  Matrix B and rank k
% Output: Factors Q and F defining rank-k approximation Bhat = Q*F',
%         set of pivots S

[m,n] = size(B);
Q = zeros(m,k);                     % Left factor (orthonormal cols)
F = zeros(n,k);                     % Right factor
S = zeros(k,1);                     % Pivots

for i = 1:k
    % Random sample using current diagonal as sampling weights
    [~,s] = datasample(1:n,1,"Weights",sqcolnorms(B));
    S(i) = s;                       % Set pivot
    Q(:,i) = B(:,s) / norm(B(:,s)); % Normalize column
    F(:,i) = B' * Q(:,i);           % Column of factor matrix
    B = B - Q(:,i) * F(:,i)';       % Modified Gram-Schmidt
end

end