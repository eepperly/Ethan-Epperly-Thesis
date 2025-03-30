function QtX = apply_Qt(hous,tau,X)
% Input:  Compact representation [hous,tau] for QR factorization, 
%         input matrix X
% Output: Product QtX = Q'*X

QtX = matlab.internal.decomposition.applyHouseholder(hous, tau,...
    X, true, size(hous,2));

end