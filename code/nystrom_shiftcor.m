function [U,D] = nystrom_shiftcor(A,n,k)
% Input:  Function A() computing matrix products A(X) = A*X, number 
%         of rows n, rank k
% Output: Shift-corrected Nystrom approximation Ahat = U*D*U', 
%         represented by orthonormal factor U and nonnegative 
%         diagonal factor D

[F,mu] = nystrom(A,n,k); % Compute (shifted) Nystrom approximation
[U,S,~] = svd(F,"econ"); % Economy-size SVD
D = max(S.^2 - mu, 0);   % Apply shift correction

end