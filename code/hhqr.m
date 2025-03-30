function [hous,tau] = hhqr(Y)
% Input:  Matrix Y
% Output: Compact representation [hous,tau] for QR factorization of Y

[hous,tau] = matlab.internal.decomposition.compactQR(Y);

end