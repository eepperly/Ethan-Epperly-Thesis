function X = spir(B,c,St,iter1,iter2)
% Input:  Matrix B, right-hand side c, sketching matrix St, and 
%         iteration count iter
% Output: Approximate least-squares solution x
    
[~,D,V] = svd(St*B,"econ");        % Sketch and compute SVD
prec = V / D;                      % Form preconditioner
y0 = Ut*(St*c);                    % Sketch-and-solve initialize
Y1 = mylsqr(@(z) B*(prec*z), ...   % Sketch-and-precondition
           @(z) prec'*(B'*z), ... 
           c, y0, iter1);
Y2 = mylsqr(@(z) B*(prec*z), ...   % Refinement step
           @(z) prec'*(B'*z), ... 
           c, Y(:,end), iter2);
X = prec * [Y1 Y2(:,2:end)];       % Change of variables

end