function X = mypcg(A,pre,b,iter)
% Input:  Function A() computing matrix products A(z) = A*z, function
%         pre applying the preconditioner, right-hand side b, and
%         number of iterations iter
% Output: Approximate solutions stacked columnwise as X

X = zeros(size(b,1),iter+1);        % Initialize solutions
z = pre(b); p = z;                  % Initialize search direction
for i = 1:iter
    v = A(p);                       % Apply A
    zb = z'*b; eta = zb / (v'*p);
    X(:,i+1) = X(:,i) + eta*p;      % Update solution
    b = b - eta*v;                  % Update residual
    z = pre(b);                     % Apply preconditioner
    gamma = z'*b/zb;
    p = z + gamma*p;                % Update search direction
end

end
