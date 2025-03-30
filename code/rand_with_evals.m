function A = rand_with_evals(evals)
% Input:  Eigenvalues evals
% Output: Hermitian matrix A with specified eigenvalues

n = length(evals);
[Q,~] = qr(randn(n));
A = Q * diag(evals) * Q';
A = (A+A')/2;

end