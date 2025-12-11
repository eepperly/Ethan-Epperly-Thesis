function [g,beta,S] = rpcholesky_active_krr(y,kernel,D,lamb,k)
% Input:  Outputs y, kernel function kernel(x,x'), inputs D, 
%         regularization lamb, and number of points k
% Output: Kernel interpolant g, coefficients beta, and subset S

Acol = @(S) kernel(D,D(S,:));            % Column subroutine
Asub = @(S) kernel(D(S,:),D(S,:));       % Submatrix subroutine
d = ones(size(D,1),1);                   % Diagonal of kernel matrix
b = ceil(k/2);                           % Block size for acc. RPC
[~,S] = acc_rpcholesky(Acol,Asub,d,k,b); % Subset selection by RPC

[g,beta] = krr(y,kernel,D(S,:),lamb);    % Apply KRR

end
