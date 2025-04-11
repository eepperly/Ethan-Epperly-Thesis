function d = bks(B,n,s)
% Input:  Function B() computing matrix products B(X) = B*X,
%         number of columns n and number of matvecs s
% Output: Estimate d of the diagonal of B

Om = random_signs(n,s);    % Matrix of random signs
d = mean(B(Om) .* Om,2);   % BKS diagonal estimate

end