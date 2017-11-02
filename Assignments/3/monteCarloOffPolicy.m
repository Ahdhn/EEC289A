function [ordinarySampling , weightedSampling] = monteCarloOffPolicy(nEpisodes)
    initialState = [true, 13, 2];
    sumOfImportanceRatio = 0;
    sumOfRewards = 0;
    
    for x = 1:nEpisodes
        [blah, reward,playerTrajectory ] = play(1, initialState);
        importanceRatioAbove = 1.0;
        importanceRatioBelow = 1.0;
        [len, blah] = size(playerTrajectory);
        
        for y=1: len
            action = playerTrajectory(y,1);
            usableAce =playerTrajectory(y,2);
            playerSum = playerTrajectory(y,3);
            dealerCard = playerTrajectory(y,4);
            
            if action == targetPolicyPlayer(playerSum)
                importanceRatioBelow = importanceRatioBelow*0.5;
            else
                importanceRatioAbove = 0;
                break;
            end            
        end
        
        importanceRatio = importanceRatioAbove / importanceRatioBelow;
        sumOfImportanceRatio(end+1) = sumOfImportanceRatio(end) + importanceRatio;
        sumOfRewards(end+1) = sumOfRewards(end) + reward*importanceRatio;       
    end
    
    sumOfImportanceRatio = sumOfImportanceRatio(2:end);
    sumOfRewards = sumOfRewards(2:end);
     
    ordinarySampling = sumOfRewards ./(1:nEpisodes);
    
    for x = 1:length(sumOfImportanceRatio)
        if sumOfImportanceRatio(x)~=0
            weightedSampling(x) = sumOfRewards(x) / sumOfImportanceRatio(x);
        else
            weightedSampling(x) =0;
        end
         
    end
    
end


function target = targetPolicyPlayer(playerSum)
    global policyPlayer;
    target = policyPlayer(playerSum);
end