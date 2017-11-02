%%%%%%%ECE 289A - An Introduction to Reinforcement Learning
%%%HW 3 by Ahmed H. Mahmoud 

%%Regenerating figure 5.4


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

%Figure 5.4
trueValue = -0.27726;
nRuns = 100;
ordinarySampling = zeros(1,10000);
weightedSampling = zeros(1,10000);
for x = 1:nRuns
    [ordinarySampling_, weightedSampling_] = monteCarloOffPolicy(10000);   
    ordinarySampling = ordinarySampling + (ordinarySampling_-trueValue).^2;
    weightedSampling = weightedSampling + (weightedSampling_-trueValue).^2;
end
ordinarySampling = ordinarySampling./nRuns;
weightedSampling = weightedSampling./nRuns;
axisX = log10(1:10000);
figure;
plot(axisX, ordinarySampling, axisX, weightedSampling);
legend('Ordinary Importance Sampling','Weighted Importance Sampling');
xlabel('Episodes (10^x)');
ylabel('Mean square error');


