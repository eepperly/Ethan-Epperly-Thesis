function d = diagprod(F,G)
% Input:  Matrices F and G of the same size
% Output: The diagonal of the product d = diag(F'*G)

d = sum(conj(F).*G,1).';

end