%%%%%%%ECE 289A - An Introduction to Reinforcement Learning
%%%HW 1 by Ahmed H. Mahmoud 

%On my laptop the code takes minute to execute and generate the plot for 
%nBandits =2000 and time =1000

%%Re-generating Figure 2.4 by converting a simplified version of the code in 
%https://github.com/ShangtongZhang/reinforcement-learning-an-introduction/blob/master/chapter02/TenArmedTestbed.py
%into matlab code 

clc;
clear;
close all;

nBandits = 2000;
time =1000;


C = 2;
eps = 0.1;
stepSize=0.1;
kArm=10;
initial_val=0;

for i=1:nBandits
    %kArm,stepSize,sampleAverages,eps, initial,UCB_param, isGradient)
    mybandUCB(i) = Bandit(kArm,stepSize,false ,0, initial_val,C, false);
    myband_eps_greedy(i) = Bandit(kArm,stepSize,false ,eps, initial_val,-1, false);
end

%combine the two instances in one vector 
[myband] = [mybandUCB; myband_eps_greedy];

%pass all instances to the bandit simulation function to get the average
%reward and (optionally) the count for choosing the best action
[bestActionCounts, averageRewards] = banditSimulation(nBandits,time,myband);

plot(1:time,averageRewards(1,:),'LineWidth',1.5);
hold on
plot(1:time,averageRewards(2,:),'LineWidth',1.5);
xlabel('Steps');
ylabel('Average Reward');
label1 = ['UCB C=' num2str(C)];
label2 = ['\epsilon - greedy, \epsilon=' num2str(eps)];
lgd = legend(label1,label2,'Location','best');
lgd.FontSize = 12;