%%%%%%%ECE 289A - An Introduction to Reinforcement Learning
%%%HW 1

%%Re-generating Figure 2.4 by converting a simplified version of the code in 
%https://github.com/ShangtongZhang/reinforcement-learning-an-introduction/blob/master/chapter02/TenArmedTestbed.py
%into matlab code 

clc;
clear;




myband = Bandit(10,0.1,false ,0.1, 0);


eps_greedy(2000,1000);
ucb(2000,1000);