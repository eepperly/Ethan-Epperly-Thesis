function s = rejection_sample(sample_proposal,pi,tau)
% Input:  Function sample_proposal() generating a sample from 
%         proposal distribution, functions pi(p) and tau(p) defining  
%         (unnormalized) density functions for proposal and target
% Output: Sample s from proposal_distribution

while true
    s = sample_proposal();     % Sample from proposal
    if rand() < tau(s) / pi(s) % Accept with probability tau(s)/pi(s)
        break
    end
end

end