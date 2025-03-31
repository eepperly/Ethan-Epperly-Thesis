function Q = haarorth(m,n)
% Input:  Dimensions m and n
% Output: Haar-random matrix Q with orthonormal columns

    [Q,R] = qr(randn(m,n),"econ"); % Orthonomalize a Gaussian matrix
    Q = Q*diag(sign(diag(R)));     % Rescale as appropriately

end