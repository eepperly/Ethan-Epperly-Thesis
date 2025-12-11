function [g,Betas] = rpcholesky_precon(y,kernel,D,lamb,k,iter)
% Input:  Outputs y, kernel function kernel(x,x'), inputs D, 
%         regularization lamb, rank k, and number of iterations iter
% Output: Kernel interpolant g and sequence Beta of coefficients
%         produced by PCG, stacked columnwise

A = kernel(D,D);                    % Generate full kernel matrix
Acol = @(S) A(:,S);                 % Column generation subroutine
Asub = @(S) A(S,S);                 % Submatrix generation subroutine
d = diag(A);                        % Diagonal of kernel matrix
b = ceil(k/2);                      % Block size for accelerated RPC
F = acc_rpcholesky(Acol,Asub,d,k,b);% LRA by accelerated RPCholesky
[U,S,~] = svd(F,"econ");            % Economy-size SVD

% Define matrix-vector product, preconditioner
matvec = @(z) A*z + lamb*z;
pre = @(z) U*((S^2 + lamb*eye(k))\(U'*z)) + (z - U*(U'*z)) / lamb;

Betas = mypcg(matvec,pre,y,iter);   % Coefficients by PCG
g = @(X) kernel(X,D) * Betas(:,end);% Define interpolant

end
