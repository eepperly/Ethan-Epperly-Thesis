function St = sparse_iid(d,m,zeta)
% Input:  Dimensions d and m, expected nnz per column zeta
% Output: Sketching matrix St

    mask = rand(d, m) < zeta/d;      % Generate mask
    St = mask .* random_signs(d,m);  % Apply mask
    St = sparse(St) / sqrt(zeta);    % Scale, convert to sparse

end