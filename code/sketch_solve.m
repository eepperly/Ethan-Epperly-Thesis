function x = sketch_solve(B,c,St)
% Input:  Matrix B, right-hand side c, and sketching matrix St
% Output: Approximate least-squares solution x

x = (St*B) \ (St*c);

end