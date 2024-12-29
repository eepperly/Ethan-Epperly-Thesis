function [U,D] = nystrom(A,n,k)
% Input:  Function A() computing matrix products A(X) = A*X, number 
%         of rows n, rank k
% Output: Shift-corrected Nystrom approximation Ahat = U*D*U', 
%         represented by orthonormal factor U and nonnegative 
%         diagonal factor D

Om = randn(n,k);                % Gaussian test matrix
Y = A(Om);                      % Matrix product Y = A*Om
mu = eps*norm(Y,"fro")/sqrt(n); % Compute shift
Y = Y + mu * Om;                % Apply shift to Y
H = Om'*Y;
R = chol((H+H')/2);             % Explicitly symmetrize H to be safe
F = Y/R;                        % Triangular substitution
[U,S,~] = svd(F,"econ");        % Economy-size SVD
D = max(S.^2 - mu, 0);          % Apply shift correction

end