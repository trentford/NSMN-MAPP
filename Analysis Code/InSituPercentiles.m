function [ obsPerc ] = InSituPercentiles(obs)
%Computes percentiles from daily VWC 
%   This function takes in daily VWC from a number of in situ stations,
%   averages across those stations, then computes percentiles from the
%   resultant dataset

    obsPerc = [];
    for i = 1:12
        sub = obs(obs(:,2) == i, :);
        nm = sub(:,5); nm(isnan(nm(:,1)) == 1, :) = [];
        [f,x] = ecdf(nm);
        
        perc(:,1:4) = sub(:,1:4);
        for ii = 1:length(sub)
            if isnan(sub(ii,5)) == 0
                [c,ind] = find(x==sub(ii,5));
                perc(ii,5) = f(c(1));
                clear c ind
            else
                perc(ii,5) = NaN;
            end
        end
        
        obsPerc = [obsPerc; perc];
        
        clear sub nm f x perc
    end
    
end

