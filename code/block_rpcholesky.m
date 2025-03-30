function [F,S] = block_rpcholesky(Acol,d,k,b)
% Input:  Function Acol for producing columns Acol(I) = A(:,I) of A,
%         diagonal d of A, rank k, block size b
% Output: Factor F defining a rank-k approximation Ahat = F*F', pivot
%         set S

n = length(d);                % Matrix size
F = zeros(n,k);               % To store output 
S = zeros(k,1);               % To store pivots
shift = 4*max(d)*eps;         % Shift to ensure positive definiteness
i = 0;                        % Index to store current position
while i < k
    b = min(b,k-i);                         % At most k pivots
    % Random sample using current diagonal as sampling weights
    Snew = datasample(1:n,b,"Weights",d,"Replace",false);
    S(i+1:i+b) = Snew;                      % Update pivots
    G = Acol(Snew) - F(:,1:i)*F(Snew,1:i)'; % Columns of residual
    R = chol(G(Snew,:) + shift*eye(b));     % Shift for stability
    F(:,i+1:i+b) = G / R;                   % Update factor
    d = d - sqrownorms(F(:,i+1:i+b));       % Update diagonal
    d = max(d,0); % Ensure nonnegative diagonal in floating point
    i = i + b;                              % Update index
end

end