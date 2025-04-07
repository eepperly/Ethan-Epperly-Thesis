function [S,R] = rejection_rpcholesky(kernel,sample_diagonal,k,d)
% Input:  Kernel function kernel(S,T) giving m*n pairwise evaluations
%         of the kernel function on the rows of an m*d matrix S and
%         n*d matrix T, function sample_diagonal() giving a sample
%         from kernel(x,x) dmu(x) outputted as a row vector,
%         number of landmarks k, dimension d of data
% Output: Set S of k selected landmarks stored as a k*d array,
%         Cholesky factor R = chol(kernel(S,S))

S = zeros(k,d);                       % Landmark set
R = zeros(k,k);                       % Cholesky factor
for i = 0:(k-1)
    while true
        S(i+1,:) = sample_diagonal(); % Proposal

        % Quantities needed to evaluate ratio
        c = R(1:i,1:i)' \ kernel(S(1:i,:),S(i+1,:));
        kss = kernel(S(i+1,:),S(i+1,:));
        d = kss - norm(c)^2;
        if rand() < d / kss       % Accept or reject
            % Update Cholesky factor upon acceptance
            R(1:i,i+1) = c;
            R(i+1,i+1) = sqrt(d);
            break                 % Exit loop and sample next point
        end
    end
end