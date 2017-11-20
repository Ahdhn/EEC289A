function temporalDifference(step, alpha)
    %n step TD method 
    global  START_STATE END_STATE GAMMA currentStateValues
    currentState = START_STATE;  
    states = currentState; %stores states for an eps
    rewards =0;%stores rewards for an eps
    time =1;
    T =inf;%length of this eps
    while true
        time = time+1;
        if time <T
            %choose action randomly
            if nbinrnd(1,0.5) == 1
                newState = currentState +1;
            else
                newState = currentState -1;
            end
            if newState == 1
                reward = -1;
            elseif  newState == 21
                reward = 1;
            else 
                reward = 0;
            end
            
            %store new state and new reward
            states(end + 1) = newState;
            rewards(end + 1)= reward;
            
            if ismember(newState, END_STATE)
                T = time;
            end
        end
        
        updateTime = time - step;
        if updateTime >=1
            returns =0;
            
            %calc corresponding rewards
            for t = updateTime : min(T,updateTime+step)  
                returns = returns + GAMMA^(t-updateTime-1)*rewards(t);
            end
            
            %add state value to the return
            if updateTime + step <= T                 
                returns = returns + GAMMA^step *currentStateValues(states(updateTime+step));
            end
            stateToUpdate = states(updateTime);
            %update the state value
            if ismember(stateToUpdate, END_STATE)==0
                currentStateValues(stateToUpdate) = currentStateValues(stateToUpdate) + alpha*(returns - currentStateValues(stateToUpdate));
            end            
        end
              
        
        if updateTime == T
            break;
        end
        currentState = newState;
    end
    
    
end