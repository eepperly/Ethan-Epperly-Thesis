function d = xnysdiag(A,n,s)
% Input:  Function A() computing matrix products A(X) = A*X, number
%         of rows n, and number of matvecs s
% Output: Estimate d of diag(A)

[F,mu,Om,R] = nystrom(A,n,s);                  % Nystrom approx 
Z = (F/R') .* (sqrownorms(inv(R)) .^ (-1/2))'; % Downdate it
d = sqrownorms(F) + (-sqrownorms(Z) + sum((Z...% Compute estimator
    .* conj(Om)) .* diagprod(Om,Z).',2))/s ...
    - mu * ones(n,1);

end