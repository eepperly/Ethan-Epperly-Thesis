function [B,c,x,r,s,V] = random_ls_problem(m,n,condB,rnorm)
% Input:  Dimensions m and n, condition number condB, and residual
%         norm rnorm
% Output: Matrix B, right-hand side c, solution x = B\c, residual r
%         = c - B*x, singular values s, and right singular vectors V

    U = haarorth(m,n+1);             % Haar-random left sing vecs
    V = haarorth(n,n);               % Haar-random right sing vecs
    s = logspace(-log10(condB),0,n); % Log-spaced sing vals
    B = U(:,1:n)*diag(s)*V';
    x = orth(randn(n,1));            % Unit-vector solution
    r = U(:,end) * rnorm;            % Residual orthog to range(B)
    c = B*x + r;

end