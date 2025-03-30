function [F,S] = nuclear_chol(A,k)
% Input:  Psd matrix A, rank k
% Output: Factor F defining a rank-k approximation Ahat = F*F', pivot
%         set S

F = zeros(size(A,1),k);                       % To store output 
S = zeros(k,1);                               % To store pivots
for i = 1:k
    nuc = sqcolnorms(A) ./ diag(A);           % Nuclear scores
    nuc(S(1:(i-1))) = 0;                      % Zero selected scores
    [~,S(i)] = max(nuc);                      % Largest nucl score 
    F(:,i) = A(:,S(i)) / sqrt(A(S(i),S(i)));  % Build factor
    A = A - F(:,i)*F(:,i)';                   % Update residual
end

end