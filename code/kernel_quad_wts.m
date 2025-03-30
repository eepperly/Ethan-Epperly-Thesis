function [w,integrator] = kernel_quad_wts(kernel,S,Au)
% Input:  Kernel function kernel(x,x'), landmark set S, and integrals
%         Au = int kernel(S,x) u(x) dx
% Output: Ideal weights w and function integrator(f) computing the
%         ideally weighted quadrature approximation w'*f(S)

w = kernel(S,S) \ Au;
integrator = @(f) w'*f(S);

end