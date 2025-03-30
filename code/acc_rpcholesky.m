function [F,S] = acc_rpcholesky(Acol,Asub,d,k,b)
% Input:  Function Acol for producing columns Acol(I) = A(:,I) of A,
%         function Asub for producing submatrices Asub(I) = A(I,I) of
%         A, diagonal d of A, rank k, block size b
% Output: Factor F defining a rank-k approximation Ahat = F*F', pivot
%         set S

n = length(d);                % Matrix size
F = zeros(n,k);               % To store output 
S = zeros(k,1);               % To store pivots
i = 0;                        % Index to store current position
while i < k
    % Random sample using current diagonal as sampling weights
    Sp = datasample(1:n,b,"Weights",d,"Replace",true);
    H = Asub(Sp) - F(Sp,1:i)*F(Sp,1:i)';    % Residual submatrix
    [T,L] = rejection_sample_submatrix(H,diag(H),k-i);
    T = Sp(T);                              % Get selected pivots
    l = length(T);                          % Number of pivots
    S(i+1:i+l) = T;                         % Update pivots
    G = Acol(T) - F(:,1:i)*F(T,1:i)';       % Columns of residual
    F(:,i+1:i+l) = G / L';                  % Update factor
    d = d - sqrownorms(F(:,i+1:i+l));       % Update diagonal
    d = max(d,0); % Ensure nonnegative diagonal in floating point
    i = i + l;                              % Update index
end

end