function [Q,F,S] = acc_rpqr(B,k,b)
% Input:  Matrix B, rank k, block size b
% Output: Factors Q and F defining rank-k approximation Bhat = Q*F',
%         set of pivots S

[m,n] = size(B);
hous = zeros(m,k); tau = zeros(k,1);% Compact Householder QR
S = zeros(k,1);                     % Pivots
i = 0;                              % Index to store current position
while i < k
    % Random sample using squared column norms as sampling weights
    Sp = datasample(1:n,b,"Weights",sqcolnorms(B(i+1:end,:)),"Replace",true);
    % Form Gram matrix of selected pivots
    Bp = B(i+1:end,Sp);           % Bottom part of B stores residual
    H = Bp'*Bp;                   % Gram matrix
    T = rejection_sample_submatrix(H,diag(H),k-i);
    T = Sp(T);                    % Get selected pivots
    l = length(T);                % Number of pivots
    S(i+1:i+l) = T;               % Update pivots
    % Update Householder QR factorization
    hous(:,i+1:i+l) = B(:,T); 
    [hous(i+1:end,i+1:i+l),tau(i+1:i+l)]... % Compute Householder
        = hhqr(hous(i+1:end,i+1:i+l));      %    reflectors
    B(i+1:end,:) = apply_Qt(hous(i+1:end,i+1:i+l),... % Update 
        tau(i+1:i+l),B(i+1:end,:));                   %   residual
    i = i + l;
end

Q = get_Q(hous,tau);
F = B(1:k,:)';

end