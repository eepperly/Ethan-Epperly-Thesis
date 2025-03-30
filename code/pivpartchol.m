function F = pivpartchol(Acol,n,s)
% Input:  Function Acol for producing columns Acol(i) = A(:,i) of A,
%         size n of A, list s = [s(1) ... s(k)] of k pivots to
%         eliminate
% Output: Factor F defining a rank-k approximation Ahat = F*F'

F = zeros(n,length(s));                       % To store output 
for i = 1:length(s)
    ai = Acol(s(i)) - F(:,1:i-1)*F(i,1:i-1)'; % ith col of A - F*F'
    F(:,i) = ai / sqrt(ai(s(i)));             % Rescale
end

end