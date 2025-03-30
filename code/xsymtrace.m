function [tr,est] = xsymtrace(A,n,s)
% Input:  Function A() computing matrix products A(X) = A*X, number
%         of rows n, and number of matvecs s
% Output: Estimate tr of trace(A), estimate est of the error
%         abs(tr - trace(A))

% Nystrom approximation
Om = random_signs(n,s);         % Test matrix of random signs
Y = A(Om);                      % Matrix product Y = A*Om
H = Om'*Y;
[L,D] = ldl((H+H')/2);          % LDLt factorization
F = Y/L';                       % Triangular substitution

% Downdating
Z = (F/D)/L;
d = diag(inv(H)); % Downdated approx is F*D*F'-Z(:,i)*Z(:,i)'/d(i)

% Compute vector of estimates
tr_vec = trace(D\(F'*F)) * ones(s,1) - sqcolnorms(Z) ./ d ...
    + abs(diagprod(Om,Z)) .^ 2 ./ d;
tr = mean(tr_vec);           % Trace estimate
est = std(tr_vec) / sqrt(s); % Error estimate

end

