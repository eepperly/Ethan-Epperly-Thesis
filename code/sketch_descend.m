function X = sketch_descend(B,c,S,alpha,beta,iter)
% Input:  Matrix B, right-hand side c, sketching matrix S, damping 
%         and momentum coefficients alpha and beta, and iteration
%         count iter
% Output: Approximate least-squares solution x
    
    X = zeros(size(B,2),iter+1);    % Initialize history
    [Q,R] = qr(S'*B,"econ");        % Sketch and compute QR       
    X(:,1) = R\(Q'*(S'*c));         % Sketch-and-solve initialize
    X(:,2) = X(:,1) + alpha * (R\(R'\(B'*(c-B*X(:,1))))); % Step one
    for i = 2:iter
        X(:,i+1) = X(:,i)...        
            + alpha * (R\(R'\(B'*(c-B*X(:,i)))))... % Damped gradient
            + beta * (X(:,i) - X(:,i-1));           % Momentum
    end

end