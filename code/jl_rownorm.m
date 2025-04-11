function srn = jl_rownorm(B,n,s)
% Input:  Function B() computing matrix products B(X) = B*X,
%         number of columns n and number of matvecs s
% Output: Estimate srn of the squared row norms of B

Om = random_signs(n,s);    % Matrix of random signs
srn = sqrownorms(B(Om))/s; % JL row norm estimate

end