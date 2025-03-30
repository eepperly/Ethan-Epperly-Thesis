function [F,S] = iid_col_nystrom(Acol,w,k)
% Input:  Function Acol for producing columns Acol(S) = A(:,S) of A,
%         sampling weights w (often uniform w = ones(n,1)), rank k
% Output: Factor F defining a rank-k approximation Ahat = F*F',
%         selected columns S

n = length(w);                                     % Get matrix size
S = datasample(1:n,k,"Replace",false,"Weights",w); % iid sample
Y = Acol(S);                                       % Get columns
mu = eps*norm(Y,"fro")/sqrt(n);                    % Compute shift
R = chol(Y(S,:) + mu*eye(k));                      % Factorize A(S,S)
F = Y/R;                                           % Get factor

end