n = 1000;
As = {rand_with_evals(linspace(1,3,n)),...
      rand_with_evals((1:n).^(-2)),...
      rand_with_evals(0.9.^(0:(n-1))),...
      rand_with_evals(0.7.^(0:(n-1))),...
      rand_with_evals([ones(ceil(n/20),1);1e-3*ones(n-ceil(n/20),1)]),...
      rand_with_evals([ones(ceil(n/20),1);1e-8*ones(n-ceil(n/20),1)])};

s_vals = 20:20:300;

for A_idx = 1:length(As)
    A = As{A_idx};
    d = diag(A);

    xnysdiags = [];
    for s = s_vals
        s
        xnysdiags(end+1) = norm(d - xnysdiag(@(X) A*X,n,s)) / norm(d);
    end
    figure()
    semilogy(s_vals,xnysdiags)
    drawnow
end