%%%%%%%ECE 289A - An Introduction to Reinforcement Learning
%%%HW 5 - Q2 by Ahmed H. Mahmoud 
%%Regenerating figure 9.10 with reference to the python code in 
%https://github.com/ShangtongZhang/reinforcement-learning-an-introduction/blob/master/chapter09/RandomWalk.py

clc;
clear;
close all;
global N_STATES END_STATES START_STATE STEP_RANGE
%number of states except for terminal one
N_STATES = 1000;
%true state values
trueStateValues = (1/1001.0).*(-1001:2:1002);
%all states
states = 2:N_STATES+1;
%start from a centeral state
START_STATE = 500;
%terminal state
END_STATES = [1, N_STATES+2];
%actions
ACTION_LEFT = -1;
ACTION_RIGHT = 1;
ACTIONS = [ACTION_LEFT, ACTION_RIGHT];
%max stride
STEP_RANGE = 100;
episodes = 5000;
numOfTilings = 50;
tileWidth = 200;
tilingOffset = 4;

%dyamic programming to find the true state values
while true
    oldTrueStateValues = trueStateValues;
    for stateid = 1:length(states)
        state = states(stateid);
        trueStateValues(state)=0;
        for action_id = 1:length(ACTIONS)
            action = ACTIONS(action_id);
            for step = 1:STEP_RANGE              
                mystep = step*action;
                newState = state + mystep;
                newState = max(min(newState, N_STATES+2),1);
                %asynchronous update for faster convergence 
                trueStateValues(state) = trueStateValues(state) +...
                    1.0/(2.0*STEP_RANGE)*trueStateValues(newState);
            end
        end
    end   
    error = sum(abs(oldTrueStateValues - trueStateValues));
    %disp(num2str(error));
    if error < 1e-2
        break;
    end    
end
%correct the state value for terminal state to 0
trueStateValues(1)=0;
trueStateValues(end)=0;




errors = zeros(2, episodes);
stateValues = zeros(1, length(states));

%an instance of TilingsValueFunction
tilingValInst = TilingsValueFunction(numOfTilings, tileWidth, tilingOffset);
%an instance of ValueFunction 
valInst = ValueFunction(floor(N_STATES/ tileWidth));

for run =1:1
    for episode = 1:episodes
        alpha = 1.0/(episode+1);
        gradientMonteCarlo(valInst, alpha);
        for n =1:length(states)        
            stateValues(n) = value(valInst,states(n));       
        end
        errors(2,episode) = errors(2,episode) + sqrt(mean((trueStateValues(2:end-1)-stateValues).^2));
    end
    for episode = 1:episodes    
        alpha = 1.0/(episode +1);
        gradientMonteCarlo(tilingValInst, alpha);
        for n =1:length(states)
            stateValues(n) = value(tilingValInst,states(n));       
        end
        errors(1,episode) = errors(1,episode) + sqrt(mean((trueStateValues(2:end-1)-stateValues).^2));
    end
    
end

plot(errors(1,:));
hold on
plot(errors(2,:));
xlabel('Episodes');
ylabel('RMSVE (one run)');
legend('Tile Coding (50 tilings)','State aggregation (one tile)');





function [newstate, reward] = TakeAction (state, action)
    %take action at state, return reward and new state
    global STEP_RANGE END_STATES
    step = randi(STEP_RANGE+1);
    step = step * action;
    newstate = state + step;
    newstate = max(min(newstate, END_STATES(2)),END_STATES(1) );
    if newstate == END_STATES(1)
        reward = -1;
    elseif newstate == END_STATES(2)
        reward =1;
    else 
        reward =0;
    end
end
function action = getAction()
    if nbinrnd(1,0.5) == 1
        action =1;
    else 
        action =-1;
    end
end
function gradientMonteCarlo(valueFunction, alpha)
    %gradient MC 
    global START_STATE END_STATES
    currentState = START_STATE;
    trajectory = currentState;
    
    %assume gamma =1, so return is just the same as the latest reward
    reward =0;
    while ~ismember(currentState, END_STATES)
        action = getAction();
        [newState, reward] = TakeAction(currentState, action);
        trajectory(end+1) = newState;
        currentState = newState;
    end
    %gradient update for each state in this trajectory
    for stateID = 1:length(trajectory)-1
        state = trajectory(stateID);
        delta = alpha *(reward - value(valueFunction, state));
        update(valueFunction, delta, state);        
    end
end




