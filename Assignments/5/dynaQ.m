function steps = dynaQ(stateActionValues, model, maze, dynaParams)
    %play an eps for dyna-q 
    currentState = maze.START_STATE;
    steps = 0;
    while ismember(START_STATE, maze.GOAL_STATES) == false
        steps = steps+1;
        %get actions 
        action = chooseAction(surrentState, stateActionValues, maze,dynaParams);
        [newState, reward] = takeAction(maze, action);
        
        %q-learning update 
        stateActionValues(currentState(1), currentState(2), action) = stateActionValues(currentState(1), currentState(2), action)+ ...
            dynaParams.alpha*(reward + dynaParams.gamma*max(stateActionValues(newState(1),newState(2),:)) -...
            stateActionValues(currentState(1), currentState(2), action));
        
        %feed the model with experince 
        feed(model,currentState, action, newState, reward);        
        
        %samples the experince from the model        
        for t = 1:dynaParams.planningSteps
            [stateSample, actionSample, newStateSample, rewardSample] = sample(model);
            stateActionValues(stateSample(1), stateSample(2), actionSample) = stateActionValues(stateSample(1), stateSample(2), actionSample)+ ...
                dynaParams.alpha*(rewardSample + dynaParams.gamma*max(stateActionValues(newStateSample(1),newStateSample(2),:))- ...
                stateActionValues(stateSamples(1), stateSample(2), actionSample));
                        
        end
        
        currentState = newState;
        
        %check if we exceed the step limit
        if steps > maze.maxSteps
            break;
        end
        
    end
end
function ac = chooseAction(state, stateActionValues, maze, dynaParams)
    if nbinrnd(1,dynaParams.epsilon) == 1
        ac = maze.actions(randi(length(maze.actions)));
    else
        [blah, ac] = max(stateActionValues(state(1),state(2),:));        
    end
end