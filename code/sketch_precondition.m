function X = sketch_precondition(B,c,S,iter)
% Input:  Matrix B, right-hand side c, sketching matrix S, and 
%         iteration count iter
% Output: Approximate least-squares solution x
    
    [U,D,V] = svd(S'*B,"econ");        % Sketch and compute SVD
    prec = V / D;                      % Form preconditioner
    y0 = U'*(S'*c);                    % Sketch-and-solve initialize
    Y = mylsqr(@(z) B*(prec*z), ...    % Preconditioned LSQR
               @(z) prec'*(B'*z), ... 
               c, y0, iter);
    X = prec * Y;                      % Change of variables

end