function [F,S] = greedy_chol(Acol,d,k)
% Input:  Function Acol for producing columns Acol(i) = A(:,i) of A,
%         diagonal d of A, rank k
% Output: Factor F defining a rank-k approximation Ahat = F*F', pivot
%         set S

F = zeros(length(d),k);                          % To store output 
S = zeros(k,1);                                  % To store pivots
for i = 1:k
    [~,S(i)] = max(d);                           % Largest diag entry 
    as = Acol(S(i)) - F(:,1:i-1)*F(S(i),1:i-1)'; % sth col of A-F*F'
    F(:,i) = as / sqrt(as(S(i)));                % Rescale
    d = d - abs(F(:,i)).^2;                      % Update diagonal
    d = max(d,0); % Ensure nonnegative diagonal in floating point
end

end