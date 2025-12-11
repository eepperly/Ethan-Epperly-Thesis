function [Q,T] = lanczos(M,z,k)
% Input:  Function M() computing matrix-vector products, vector z,
%         and number of steps k
% Output: Matrix Q with orthonormal columns and tridiagonal matrix T

Q = zeros(length(z),k+1); T = zeros(k+1);    % Initialize
Q(:,1) = z / norm(z);                        % Normalize first vector
for i = 1:k
    if i == 1
        Q(:,2) = M(Q(:,1));                       % apply matrix
    else
        Q(:,i+1) = M(Q(:,i)) - T(i,i-1)*Q(:,i-1); % 3-term recurrence
    end
    T(i,i) = Q(:,i+1)'*Q(:,i);
    Q(:,i+1) = Q(:,i+1) - T(i,i)*Q(:,i);          % orthogonalize
    T(i,i+1) = norm(Q(:,i+1)); T(i+1,i) = T(i,i+1); 
    Q(:,i+1) = Q(:,i+1) / T(i,i+1);               % normalize
end
T = T(1:k,1:k); Q = Q(:,1:k);                % Shrink to size k

end