function [bestActionCounts, averageRewards] = banditSimulation (nBandits, time, bandits)
    %simluate the k-arm bandits game using the parameters specified in
    %array bandits (member of Bandit class) number of times = time with
    %nBandits steps 
        
    %initialize the best actions and average rewards 
    [bandits_Num, blah] = size(bandits);
    bestActionCounts = zeros(bandits_Num,time);
    averageRewards = zeros(bandits_Num,time);
    
    for banditID = 1:bandits_Num %for every configuration passed 
        for i=1:nBandits %for every step
            for t = 1:time %to average over this step
                action = getAction(bandits(banditID,i));
                reward = takeAction(bandits(banditID,i),action);
                
                %accumelate for (later) calc of average rewards at this
                %time 
                averageRewards(banditID,t) = averageRewards(banditID,t) + reward; 
                
                %increment if we choose the best action
                if action == bandits(banditID,i).bestAction
                    bestActionCounts(banditID,t) = bestActionCounts(banditID,t)+1;
                end
            end 
        end
        %don't know why it doesn't work without the colon notation 
        bestActionCounts(banditID,:) = bestActionCounts(banditID,:)/nBandits;
        averageRewards(banditID,:) = averageRewards(banditID,:)/nBandits;        
    end
end 