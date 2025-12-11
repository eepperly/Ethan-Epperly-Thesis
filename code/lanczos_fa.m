function fMZ = lanczos_fa(M,Z,k,f)
% Input:  Function M() computing matrix-vector products, matrix Z,
%         number of steps k, function f
% Output: Approximation fMZ to f(M)*Z

fMZ = zeros(size(Z));                                % initialize
for i = 1:size(Z,2)
    [Q,T] = lanczos(M,Z(:,i),k);                     % run Lanczos
    [V,d] = eig(T,"vector");                         % diagonalize
    fMZ(:,i) = norm(Z(:,i))*(Q*(V*(f(d).*V(1,:)'))); % Lanczos-FA
end

end