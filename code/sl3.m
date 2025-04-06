function srn = sl3(B,Bt,m,n,s)
% Input:  Functions B() and Bt() computing matrix products B(X) = 
%         B*X and Bt(X) = B'*X, dimensions m and n, and number of
%         matvecs s
% Output: Estimate srn of the squared row norms of B

k = floor(s/3);                             % Approximation rank s/3
Om = random_signs(m,k);                     % Matrix of random signs
Y = Bt(Om);                                 % Matvecs with B'
[Q,~] = qr(Y,"econ");                       % Orthonormalize                         
BQ = B(Q);                                  % Matvecs with B
k = s - 2*k;                                % Remaining matvecs
Ga = random_signs(n,k);                     % Matrix of random signs
Ga = Ga - Q * (Q'*Ga);                      % Orthogonalize against Q
srn = sqrownorms(BQ) + sqrownorms(B(Ga))/k; % SL3 estimator

end