function [T,L] = rejection_sample_submatrix(H,u,lmax)
% Input:  Psd matrix H and proposal distribution u
% Output: Indices T of selected pivots, Cholesky factor L
%         = chol(H(T,T), "lower")

b = size(H,1);                         % Matrix dimension
T = zeros(b,1);                        % Buffer for accepted pivots
num_accepts = 0;                       % Counter for acceptances

for j = 1:b
    if u(j) * rand() > H(j,j)            % Rejection sampling
        continue                         % If reject, continue
    end
    num_accepts = num_accepts + 1;       % Increment counter
    T(num_accepts) = j;                  % Accept i as pivot
    H(j:b,j) = H(j:b,j) / sqrt(H(j,j));  % Overwrite H with Cholesky
    % Update residual (Schur complement)
    H(j+1:b,j+1:b) = H(j+1:b,j+1:b) - H(j+1:b,j)*H(j+1:b,j)';
    if num_accepts == lmax; break; end   % At most lmax pivots
end

T = T(1:num_accepts);                    % Fix buffer size
L = tril(H(T,T));                        % Extract Cholesky factor

end