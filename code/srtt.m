classdef srtt
    properties
        signs  % Matrix of random signs
        idx    % Indexes for subsampling
    end
    
    methods
        function St = srtt(d, m, rounds)
            St.signs = random_signs(m, rounds);  % Generate signs 
            St.idx = randsample(m, d, false);    % Subsample indices
        end
        
        function y = mtimes(St, x)
            for i = 1:size(St.signs,2)
                x = dct(St.signs(:,i) .* x);     % Random trig trans
            end
            y = x(St.idx,:);                     % Subsample
            y = sqrt(length(x)/length(y)) * y;   % Rescale
        end
    end
end