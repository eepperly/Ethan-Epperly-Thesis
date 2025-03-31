function be = backerr_est(B,c,x,s,V)
% Input:  Matrix B, column vector of singular values s,
%         right singular vectors V of B, right-hand side c, and
%         approximate LS solution x
% Output: Backward error estimate be for min ||c-B*x||
    
    r = c - B*x;                  % Residual
    om = norm(r) / norm(x);
    be = om / norm(r) * norm((V'*(B'*r)) ./ (s.^2 + om^2).^0.5);

end