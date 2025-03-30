function [S,T,U,G,P] = importance_cur(B,k,l)
% Input:  Matrix B and ranks k and l
% Output: Column and row sets S and T, k*l core matrix U defining CUR
%         approximation B ~ B(:,S) * U * B(T,:), matrices G and P 
%         defining (more) stable representation U = P \ G

[m,n] = size(B);
S = datasample(1:n,k,"Replace",false,"Weights",sqcolnorms(B));
T = datasample(1:m,l,"Replace",false,"Weights",sqrownorms(B));
[Q,R] = qr(B(T,:)',"econ");
W = (B*Q)/R';                          % Row interpolation matrix
[Q,P] = qr(B(:,S),"econ");
G = Q'*W;                              % Well-conditioned matrix
U = P \ G;                             % Standard core matrix U

end