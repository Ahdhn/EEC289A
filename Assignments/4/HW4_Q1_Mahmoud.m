%%%%%%%ECE 289A - An Introduction to Reinforcement Learning
%%%HW 4 - Q1 by Ahmed H. Mahmoud 
%%Regenerating figure 6.5 
clc;
clear;
close all;

global actionRewards startState goalState EPSILON actions GAMMA ...
       actionDestination stateActionValuesSarsa stateActionValuesQLearning

WORLD_HEIGHT = 4; %world height
WORLD_WIDTH = 12; %world width
EPSILON = 0.1; %probability for exploration
ALPHA = 0.5; %step size
GAMMA = 1; %gamma for Q-Learning and Expected Sarsa
ACTION_UP = 1; %all possible actions
ACTION_DOWN = 2;
ACTION_LEFT = 3;
ACTION_RIGHT = 4;
actions = [ACTION_UP, ACTION_DOWN, ACTION_LEFT, ACTION_RIGHT];

%initial state action pair values
stateActionValues = zeros(WORLD_HEIGHT, WORLD_WIDTH, 4);
startState = [4, 1];
goalState = [4, 12];

%reward for each action in each state
actionRewards = zeros(WORLD_HEIGHT, WORLD_WIDTH, 4);
actionRewards(:, :, :) = -1.0;
actionRewards(3, 2:11, ACTION_DOWN) = -100.0;
actionRewards(4, 1, ACTION_RIGHT) = -100.0;

%set up destinations for each action in each state
for x =1:WORLD_HEIGHT
    for y = 1:WORLD_WIDTH
        actionDestination(x,y).ACTION_UP =    [max(x-1,1), y];
        actionDestination(x,y).ACTION_LEFT =  [x max(y-1,1)];
        actionDestination(x,y).ACTION_RIGHT = [x, min(y+1, WORLD_WIDTH)];
        if x == 3 && 2<=y && y<=11
            actionDestination(x,y).ACTION_DOWN = startState;
        else
            actionDestination(x,y).ACTION_DOWN = [min(x + 1, WORLD_HEIGHT), y];
        end        
    end
end
actionDestination(4,1).ACTION_RIGHT = startState;



averageRange = 10;% averaging the reward sums from 10 successive episodes
runs = 20;%perform 20 independent runs
nEpisodes = 500;%episodes of each run

rewardsSarsa = zeros(1,nEpisodes);
rewardsQLearning = zeros(1,nEpisodes);

for x = 1:runs
    stateActionValuesSarsa = stateActionValues;
    stateActionValuesQLearning = stateActionValues;
    for y=1:nEpisodes
        %cut off the value by -100 to draw the figure more elegantly
        rewardsSarsa(y) = rewardsSarsa(y) + max(sarsa(@chooseAction,@actionDestinationMapping,ALPHA), -100); 
        rewardsQLearning(y) = rewardsQLearning(y) + max(qLearning(@chooseAction,@actionDestinationMapping, ALPHA), -100);
    end   
end
%averaging over runs
rewardsSarsa = rewardsSarsa/runs;
rewardsQLearning = rewardsQLearning / runs;
    
%avergaing over successive episodes 
smoothedRewardsSarsa = smooth(rewardsSarsa,10);
smoothedRewardsQLearning =smooth(rewardsQLearning,10);


plot (1:nEpisodes, smoothedRewardsSarsa,1:nEpisodes,smoothedRewardsQLearning);
legend('Sarsa','Q-Learning','Location','southeast');
xlabel('Episodes','FontSize', 13,'FontWeight','bold');
ylabel('Sum of rewards during episodes','FontSize', 13,'FontWeight','bold');


function action = chooseAction(state, stateActionValues)
    %choose an action based on epsilon greedy algorithm
    global EPSILON actions
    if nbinrnd(1,EPSILON) == 1
        action = actions(randi(length(actions)));
    else
        [blah, action] = max(stateActionValues(state(1),state(2),:));
    end
    if action <1 || action >4
        error('noooo');
    end
end

function state = actionDestinationMapping(x,y,Action)
    %map Action (numeric) to one of the actionDestination (char)
    %return the state 
    global actionDestination;
    if Action == 1 %ACTION_UP
        state = actionDestination(x,y).ACTION_UP;
    elseif Action ==2 %ACTION_DOWN
        state = actionDestination(x,y).ACTION_DOWN;
    elseif Action ==3%ACTION_LEFT
        state = actionDestination(x,y).ACTION_LEFT;
    elseif Action ==4%ACTION_RIGHT
        state = actionDestination(x,y).ACTION_RIGHT;
    else
        error('Wrong index to actionDestination!!!');
    end    
end
