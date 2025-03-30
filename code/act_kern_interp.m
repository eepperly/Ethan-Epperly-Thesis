function [fhat,beta,S] = act_kern_interp(f,kernel,sample_diag,k,d)
% Input:  Function f(S) giving evaluation m*1 of function on the rows
%         of m*d matrix S, kernel function kernel(S,T) giving m*n
%         pairwise evaluations of the kernel function on the rows of 
%         an m*d matrix S and n*d matrix T, function
%         sample_diag() giving a sample from kernel(x,x) dmu(x)
%         outputted as a row vector, number of interpolating points
%         k, dimension d of data
% Output: Regression function fhat(x) defined by coefficients beta,
%         set of landmarks S

% Sample landmarks
S = rpcholesky_rejection(kernel,sample_diag,k,d); 

% Get interpolation coefficients, define interpolant
beta = kernel(S,S) \ f(S);
fhat = @(X) kernel(X,S) * beta;

end