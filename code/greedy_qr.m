function [Q,F,S] = greedy_qr(B,k)
% Input:  Matrix B, rank k
% Output: Factors Q and F defining rank-k approximation Bhat = Q*F',
%         set of pivots S

m = size(B,1);
hous = zeros(m,k); tau = zeros(k,1);% Compact Householder QR
S = zeros(k,1);                     % Pivots
for i = 1:k
    % Select column with largest residual norm
    [~,S(i)] = max(sqcolnorms(B(i+1:end,:)));
    % Update Householder QR factorization
    hous(:,i) = B(:,S(i)); 
    [hous(i:end,i),tau(i)] = hhqr(hous(i:end,i));
    % Update residual
    B(i:end,:) = apply_Qt(hous(i:end,i),tau(i),B(i:end,:));
end

Q = get_Q(hous,tau);
F = B(1:k,:)';

end