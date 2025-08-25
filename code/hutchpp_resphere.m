function tr = hutchpp_resphere(B,n,s)
% Input:  Function B() computing matrix products B(X) = B*X, number
%         of rows n, and number of matvecs s
% Output: Estimate tr of trace(B)

k = floor(s/3);                           % Rank is s/3
Om = randn(n,k);                          % Random Gaussian matrix
Y = B(Om);                                % Collect matvecs           
[Q,~] = qr(Y,"econ");                     % Randomized SVD
BQ = B(Q);                                % Collect matvecs

k2 = s-2*k;                               % Remaining matvecs
Ga = randn(n,k2);                         % Gaussian random vectors
X = Ga - Q*(Q'*Ga);                       % Orthogonalize against Q
X = sqrt(n-k) * X ./ sqcolnorms(X)'.^0.5; % Resphere
BX = B(X);                                % Collect matvecs
tr = trace(Q'*BQ) + trace(X'*BX) / k2;    % Hutch++ estimator

end