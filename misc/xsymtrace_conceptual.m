function [tr,est] = xsymtrace_conceptual(A,n,s)
% Input:  Function A() computing matrix products A(X) = A*X, number
%         of rows n, and number of matvecs s
% Output: Estimate tr of trace(A), estimate est of the error
%         abs(tr - trace(A))

% Nystrom approximation
Om = random_signs(n,s);         % Test matrix of random signs
Y = A(Om);                      % Matrix product Y = A*Om
H = Om'*Y;

tr_vec = zeros(s,1);

for i = 1:s
    Yi = Y(:,[1:(i-1) (i+1):s]);
    Hi = H([1:(i-1) (i+1):s],[1:(i-1) (i+1):s]);
    Ahat = (Yi/Hi)*Yi';
    tr_vec(i) = trace(Ahat) + Om(:,i)'*(A(Om(:,i)) - Ahat*Om(:,i));
end

tr = mean(tr_vec);           % Trace estimate

end