function srn = sl4(B,Bt,m,n,s)
% Input:  Functions B() and Bt() computing matrix products B(X) = 
%         B*X and Bt(X) = B'*X, dimensions m and n, and number of
%         matvecs s
% Output: Estimate srn of the squared row norms of B

k = floor(s/4);                             % Approximation rank s/4
Om = random_signs(m,k);                     % Matrix of random signs
Y = Bt(B(Om));                              % Matvecs with B and B'
[Q,~] = qr(Y,"econ");                       % Orthonormalize    
BQ = B(Q);                                  % Matvecs with B
k = s - 3*k;                                % Remaining matvecs
Ga = random_signs(n,k);                     % Matrix of random signs
Ga = Ga - Q * (Q'*Ga);                      % Orthogonalize against Q
srn = sqrownorms(BQ) + sqrownorms(B(Ga))/k; % SL4 estimator

end