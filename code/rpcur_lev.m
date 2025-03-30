function [S,T,U,G,P] = rpcur_lev(B,k,l)
% Input:  Matrix B and ranks k and l
% Output: Column and row sets S and T, k*l core matrix U defining CUR
%         approximation B ~ B(:,S) * U * B(T,:), matrices G and P 
%         defining (more) stable representation U = P \ G

[Q,~,S] = acc_rpqr(B,k,floor(k/2)); % RPQR approximation B ~ Q*F'
lev = sqrownorms(Q);                % Leverage scores of Q
T = datasample(1:size(B,1),l,"Weights",lev); % Leverage score sample
w = lev(T) .^ (-1/2);               % Reweight rows
[Q,P,p] = qr(w .* B(T,S),"econ","vector");   % Pivoted QR
dP = abs(diag(P));                  % Diagonal of triangular factor
idx = find(dP > 20*eps*max(dP));    % Find large diagonal entries
Q = Q(:,idx); P = P(idx,idx);       % Delete neglible entries of Q,P
S = S(p(idx));                      % Filter pivots
G = Q' .* w.';                      % Well-conditioned matrix
U = P \ G;                          % Standard core matrix U

end