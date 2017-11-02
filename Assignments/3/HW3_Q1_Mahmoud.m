%%%%%%%ECE 289A - An Introduction to Reinforcement Learning
%%%HW 3 by Ahmed H. Mahmoud 

%%Regenerating figure 5.1 


clc;
clear;
close all;

global policyPlayer policyDealer;

ACTION_HIT = 0; 
ACTION_STAND = 1;
actions = [ACTION_HIT, ACTION_STAND];

%policy for player
policyPlayer = zeros(1,22);
policyPlayer(21) = ACTION_STAND;
policyPlayer(22) = ACTION_STAND;

%policy for dealer 
policyDealer = zeros(1,22);
policyDealer(18:22) = 1;

%Figure 5.1
[statesUsableAce2, statesNoUsableAce2] = monteCarloOnPolicy(10000);
plotData(statesUsableAce2, 'Usable Ace after 10K Episodes');
plotData(statesNoUsableAce2, 'No Usable Ace after 10K Episodes');

[statesUsableAce1, statesNoUsableAce1] = monteCarloOnPolicy(500000);
plotData(statesUsableAce1,'Usable Ace after 500K Episodes');
plotData(statesNoUsableAce1,'No Usable Ace after 500K Episodes');


