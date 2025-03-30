function [tr,est] = xnystrace(A,n,s)
% Input:  Function A() computing matrix products A(X) = A*X, number
%         of rows n, and number of matvecs s
% Output: Estimate tr of trace(A), estimate est of the error
%         abs(tr - trace(A))

[F,mu,Om,R] = nystrom(A,n,s);                  % Nystrom approx 
Z = (F/R') .* (sqrownorms(inv(R)) .^ (-1/2))'; % Downdate it

% Compute vector of estimates
tr_vec = norm(F,"fro")^2 * ones(s,1) - sqcolnorms(Z) ...
    + abs(diagprod(Z,Om)) .^ 2 - mu * n * ones(s,1);
tr = mean(tr_vec);           % Trace estimate
est = std(tr_vec) / sqrt(s); % Error estimate

end

