function [Q,F,S] = acc_rpqr_bgs(B,k,b)
% Input:  Matrix B, rank k, block size b
% Output: Factors Q and F defining rank-k approximation Bhat = Q*F',
%         set of pivots S

[m,n] = size(B);
Q = zeros(m,k);                     % Left (orthonormal) factor
F = zeros(n,k);                     % Right factor
S = zeros(k,1);                     % Pivots
i = 0;                              % Index to store current position
while i < k
    % Random sample using squared column norms as sampling weights
    Sp = datasample(1:n,b,"Weights",sqcolnorms(B),"Replace",true);
    H = B(:,Sp)'*B(:,Sp);                  % Gram matrix of residual
    T = rejection_sample_submatrix(H,diag(H),k-i);
    T = Sp(T);                             % Get selected pivots
    l = length(T);                         % Number of pivots
    S(i+1:i+l) = T;                        % Update pivots
    [Q(:,i+1:i+l),~] = qr(B(:,T),"econ");  % Orthonormalize selection
    F(:,i+1:i+l) = B'*Q(:,i+1:i+l);        % Update factor
    B = B - Q(:,i+1:i+l)*F(:,i+1:i+l)';    % Update residual
    B = B - Q(:,i+1:i+l)*(Q(:,i+1:i+l)'*B);% Orthogonalize twice
    i = i + l;                             % Update index
end

end