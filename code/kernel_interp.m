function [g,beta] = kernel_interp(y,kernel,D)
% Input:  Outputs y, kernel function kernel(x,x'), and inputs D
% Output: Kernel interpolant g and coefficients beta

beta = kernel(D,D) \ y;      % Get interpolation coefficients
g = @(X) kernel(X,D) * beta; % Define interpolant

end