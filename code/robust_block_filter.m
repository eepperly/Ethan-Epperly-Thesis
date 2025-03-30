function [T,L] = robust_block_filter(H,tau,lmax)
% Input:  Psd matrix H, tolerance tau, maximum number of pivots to
%         accept bmax
% Output: Set of pivots T, Cholesky factor L = chol(H(T,T), "lower")

b = size(H,1);         % Matrix dimension
F = zeros(size(H));    % Store Cholesky factor 
T = zeros(b,1);        % Store pivots
d = diag(H);           % Diagonal of H
orig_trace = sum(d);   % Original trace of H

for i = 1:min(b,lmax)  % Don't exceed lmax pivots
    [~,T(i)] = max(d);                           % Largest diag entry 
    hi = H(:,T(i)) - F(:,1:i-1)*F(T(i),1:i-1)';  % ith col of H-F*F'
    F(:,i) = hi / sqrt(hi(T(i)));                % Rescale
    d = d - abs(F(:,i)).^2;                      % Update diagonal
    if sum(d) <= tau * orig_trace; break; end    % Terminate?
end

T = T(1:i);            % Extract pivots
L = F(T,1:length(T));  % Extract Cholesky factor

end