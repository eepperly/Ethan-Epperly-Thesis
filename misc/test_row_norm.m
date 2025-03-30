n = 1000;
As = {rand_with_evals(linspace(1,3,n)),...
      rand_with_evals((1:n).^(-2)),...
      rand_with_evals(0.9.^(0:(n-1))),...
      rand_with_evals(0.7.^(0:(n-1))),...
      rand_with_evals([ones(ceil(n/20),1);1e-3*ones(n-ceil(n/20),1)]),...
      rand_with_evals([ones(ceil(n/20),1);1e-8*ones(n-ceil(n/20),1)])};

s_vals = 10:10:300;

for A_idx = 1:length(As)
    A = As{A_idx};
    srn = vecnorm(A,2,2).^2;

    xrows = [];
    xsymrows = [];
    sl3s = [];
    sl4s = [];
    for s = s_vals
        s
        xrows(end+1) = norm(srn - xrownorm(@(X) A*X, @(X) A*X, n, s)) / norm(srn);
        xsymrows(end+1) = norm(srn - xsymrownorm(@(X) A*X, n, s)) / norm(srn);
        sl3s(end+1) = norm(srn - sl3(@(X) A*X, @(X) A*X, n, n, s)) / norm(srn);
        sl4s(end+1) = norm(srn - sl4(@(X) A*X, @(X) A*X, n, n, s)) / norm(srn);
    end
    figure()
    semilogy(s_vals,xrows); hold on
    semilogy(s_vals,xsymrows);
    semilogy(s_vals,sl3s)
    semilogy(s_vals,sl4s)
    legend('XRowNorm','XSymRowNorm','SL3','SL4')
    drawnow
end