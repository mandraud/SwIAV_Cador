function [step, j, k] = GetEventSows( Rate )
% Calculate the stepsize and sample a random event based on the probabilities from 'Rate'.


    SumRates=sum(Rate(:));
    
    if SumRates>0,
        step = -log(rand(1))/SumRates;

        % find which event to do 
        m=find(cumsum(Rate(:))>=rand(1)*SumRates, 1);
        [d1, d2] = size(Rate);
        k = find(m <= (d1.*(1:d2)),1);
        j = rem(m,d1);
        j(j==0) = d1; 
    else
        %%
        step = 0;
        j = 0;
        k = 0;

    end
    
end

