function [fhat,beta] = krr(y,kernel,D,lamb)
% Input:  Data y, kernel function kernel(x,x'), data D, and 
%         regularization lamb >= 0
% Output: Kernel interpolant fhat and coefficients beta

A = kernel(D,D) + lamb*eye(size(D,1)); % kappa(D,D) + lamb*I
beta = A \ y;                          % Get coefficients
fhat = @(X) kernel(X,D) * beta;        % Define regression function

end