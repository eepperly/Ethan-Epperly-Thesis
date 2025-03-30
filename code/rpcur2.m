function [S,T,U,G,P] = rpcur2(B,k,l)
% Input:  Matrix B and ranks k and l
% Output: Column and row sets S and T, k*l core matrix U defining CUR
%         approximation B ~ B(:,S) * U * B(T,:), matrices G and P 
%         defining (more) stable representation U = P \ G

[~,F,T] = acc_rpqr(B',l,floor(l/2));  % RPQR approximation B ~ F*Q'
W = F/tril(F(T,:));                   % Row interpolation matrix
[Q,F,S] = acc_rpqr(B,k,floor(k/2));   % RPQR approximation B ~ Q*F'
P = tril(F(S,:))';                    % R factor for qr(B(:,S))
G = Q'*W;                             % Well-conditioned matrix
U = P \ G;                            % Standard core matrix U

end