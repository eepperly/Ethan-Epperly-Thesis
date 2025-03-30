function W = ppqr2id(F,S)
% Input:  Factor F defining a column projection approximation
%         Bhat = Q*F' and pivots S
% Output: Interpolation matrix W defining interpolative decomposition
%         Bhat = B(:,S)*W'

W = F / F(S,:);

end