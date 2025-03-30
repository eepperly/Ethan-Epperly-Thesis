function Q = get_Q(hous,tau)
% Input:  Compact representation [hous,tau] for QR factorization
% Output: Thin QR factor Q

[m,k] = size(hous);
Q = matlab.internal.decomposition.applyHouseholder(hous, tau,...
    [eye(k);zeros(m-k,k)], false, k);

end