function X = mylsqr(B,Bt,c,x0,iter)
% Input:  Functions B() and Bt() computing matrix products B(z) = 
%         B*z and Bt(z) = B'*z, right-hand side c, initialization x0,
%         number of iterations iter
% Output: Approximate least-squares solutions stacked columnwise as X
    
% First step of Golub-Kahan bidiagonalization
r = c - B(x0); beta = norm(r); u = r / beta;
v = Bt(u); alpha = norm(v); v = v / alpha;
w = v; phibar = beta; rhobar = alpha;

X = zeros(size(v,1),iter+1); X(:,1) = x0; % Initialize history
for i = 1:iter
    % Golub-Kahan bidiagonalization
    u = B(v) - alpha*u; beta = norm(u); u = u / beta;
    v = Bt(u) - beta*v; alpha = norm(v); v = v / alpha;
    rho = sqrt(rhobar^2 + beta^2);
    c = rhobar / rho;
    s = beta / rho;
    theta = s * alpha;
    rhobar = - c * alpha;
    phi = c * phibar;
    phibar = s * phibar;

    % Update solution, search direction
    X(:,i+1) = X(:,i) + (phi/rho) * w;
    w = v - (theta/rho) * w;
end

end