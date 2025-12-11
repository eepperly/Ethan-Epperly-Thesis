function [U,D,jack] = nystrom_jack(A,n,k,f)
% Input:  Function A() computing matrix products A(X) = A*X, number 
%         of rows n, rank k, spectral transformation f()
% Output: Shift-corrected Nystrom approximation Ahat = U*D*U', 
%         represented by factors U and D, and jackknife standard
%         deviation estimate jack

[F,mu,~,R] = nystrom(A,n,k); % Compute Nystrom approximation
[U,S,V] = svd(F,"econ");     % Economy-size SVD
D = max(S.^2 - mu, 0);       % Apply shift correction

W = ((S*V')/R')*diag(sqrownorms(inv(R)).^(-0.5)); % Downdating matrix
Xs = zeros(k,k,k);                      % Jackknife replicates
for j = 1:k
    [Q,d] = eig(D - W(:,j)*W(:,j)',"vector");
    d = d(end:-1:1); Q = Q(:,end:-1:1); % Evals in decreasing order
    Xs(:,:,j) = Q * diag(f(d)) * Q';    % Spectral transform
end
jack = norm(Xs - mean(Xs,3),"fro");     % Jackknife estimate

end