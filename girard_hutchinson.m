function tr = girard_hutchinson(B,n,s)
% Input:  Function B() computing matrix products B(X) = B*X, number
%         of rows n, and number of matvecs s
% Output: Estimate tr of trace(B)

tr = 0;
for i = 1:s
    om = random_signs(n,1);  % Generate vector of random signs
    tr = tr + om'*B(om) / s; % Update trace estimate
end

end