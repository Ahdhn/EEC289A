function rewards = sarsa(chooseAction,actionDestinationMapping, stepSize)
    %an episode with Sarsa
    
    global startState goalState actionRewards GAMMA stateActionValuesSarsa 
    currentState = startState;
    currentAction = chooseAction(currentState, stateActionValuesSarsa);
    rewards =0;
    
    while currentState(1)~= goalState(1) || currentState(2)~= goalState(2)        
        newState = actionDestinationMapping(currentState(1),currentState(2),currentAction);
        newAction = chooseAction(newState,stateActionValuesSarsa);
        reward = actionRewards(currentState(1), currentState(2), currentAction);
        rewards = rewards + reward;
        valueTarget = GAMMA*stateActionValuesSarsa(newState(1), newState(2), newAction);
        
        
        %sarsa update 
        stateActionValuesSarsa(currentState(1), currentState(2), currentAction) =...
            stateActionValuesSarsa(currentState(1), currentState(2), currentAction) + ...
            stepSize* (reward + valueTarget - stateActionValuesSarsa(currentState(1), currentState(2), currentAction));
        currentState = newState;
        currentAction = newAction;
    end
end 



