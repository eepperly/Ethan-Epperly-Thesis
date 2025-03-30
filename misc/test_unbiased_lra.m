n = 1000;
As = {rand_with_evals(linspace(1,3,n)),...
      rand_with_evals((1:n).^(-2)),...
      rand_with_evals(0.9.^(0:(n-1))),...
      rand_with_evals(0.7.^(0:(n-1))),...
      rand_with_evals([ones(ceil(n/20),1);1e-3*ones(n-ceil(n/20),1)]),...
      rand_with_evals([ones(ceil(n/20),1);1e-8*ones(n-ceil(n/20),1)])};
% As = {rand_with_evals(linspace(1,3,n)),...
%       rand_with_evals((1:n).^(-2)),...
%       rand_with_evals(0.9.^(0:(n-1))),...
%       rand_with_evals(0.7.^(0:(n-1)))};

k_vals = 10:10:150;
trials = 10;

sineangle = @(U,V) norm(V - U * (U'*V), 'fro');

close all
for A_idx = 1:length(As)
    A = As{A_idx};

    % Loss 1: Frobenius norm
    Afro = norm(A,'fro');
    loss = @(UU,SS,VV) norm(A - UU*SS*VV','fro') / norm(A,'fro');

    % % Loss 2: Left singular vector
    % [UA,SA,VA] = svd(A,"econ");
    % loss = @(UU,SS,VV) sineangle(UU(:,1:min(50,end)),UA(:,1:min(50,end)));

    % % Loss 3: Right singular vector
    % [UA,SA,VA] = svd(A,"econ");
    % loss = @(UU,SS,VV) sineangle(VV(:,1:min(50,end)),VA(:,1:min(50,end)));

    % % Loss 4: Singular value
    % [UA,SA,VA] = svd(A,"econ");
    % loss = @(UU,SS,VV) abs(SA(1,1) - SS(1,1)) / SA(1,1);

    % % Loss 5: Linear functional
    % loss = @(UU,SS,VV) abs(sum(A,"all") - sum(UU*SS*VV',"all")) / abs(sum(A,"all"));

    rsvds = [];
    usvds = [];
    for s = k_vals
        s
        rsvds(end+1) = 0;
        usvds(end+1) = 0;
        for i = 1:trials
            [U,S,V] = rsvd(@(X) A*X, @(X) A*X, 1000, floor(s/2));
            rsvds(end) = rsvds(end) + loss(U,S,V)/trials;
            [U,S,V] = usvd(@(X) A*X, @(X) A*X, 1000, floor(s/2));
            usvds(end) = usvds(end) + loss(U,S,V)/trials;
        end
    end
    figure()
    semilogy(k_vals,rsvds); hold on
    semilogy(k_vals,usvds)
    drawnow
end