function rewards = qLearning(chooseAction, actionDestinationMapping, stepSize)
    %episode with q-learning 
    global startState goalState actionRewards GAMMA stateActionValuesQLearning
    currentState = startState;
    rewards =0.0;
    while currentState(1)~= goalState(1) || currentState(2)~= goalState(2)
        currentAction = chooseAction(currentState, stateActionValuesQLearning);
        reward = actionRewards(currentState(1), currentState(2), currentAction);
        rewards = rewards + reward;
        newState = actionDestinationMapping(currentState(1), currentState(2), currentAction);
        
        %q learning update 
        r1 = stateActionValuesQLearning(currentState(1),currentState(2),currentAction);
        r2 = max(stateActionValuesQLearning(newState(1), newState(2),:));
        stateActionValuesQLearning(currentState(1), currentState(2), currentAction)=...
            stateActionValuesQLearning(currentState(1), currentState(2), currentAction)+ stepSize*( reward + GAMMA * r2- r1 );
        currentState = newState;
    end
end