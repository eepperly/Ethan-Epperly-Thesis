function X = sketch_precondition_zero(B,c,St,iter)
% Input:  Matrix B, right-hand side c, sketching matrix St, and 
%         iteration count iter
% Output: Approximate least-squares solution x
    
[~,D,V] = svd(St*B,"econ");        % Sketch and compute SVD
prec = V / D;                      % Form preconditioner
y0 = zeros(size(B,2),1);           % **Zero** initialization
Y = mylsqr(@(z) B*(prec*z), ...    % Preconditioned LSQR
           @(z) prec'*(B'*z), ... 
           c, y0, iter);
X = prec * Y;                      % Change of variables

end