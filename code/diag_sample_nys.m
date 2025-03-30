function [F,S] = diag_sample_nys(Acol,d,k)
% Input:  Function Acol for producing columns Acol(i) = A(:,i) of A,
%         diagonal d of A, rank k
% Output: Factor F defining a rank-k approximation Ahat = F*F', pivot
%         set S

S = datasample(1:length(d),k,"Replace",false,"Weights",d);
AS = Acol(S);             % Columns of A
F = AS / chol(AS(S,:));   % Factor matrix

end