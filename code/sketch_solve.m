function x = sketch_solve(B,c,S)
% Input:  Matrix B, right-hand side c, and sketching matrix S
% Output: Approximate least-squares solution x

x = (S'*B) \ (S'*c);

end