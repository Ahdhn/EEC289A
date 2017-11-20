%%%%%%%ECE 289A - An Introduction to Reinforcement Learning
%%%HW 5 - Q1 by Ahmed H. Mahmoud 
%%Regenerating figure 8.5 with reference to the python code in 
%https://github.com/ShangtongZhang/reinforcement-learning-an-introduction/blob/master/chapter08/VariousMaze.py

clc;
clear;
close all;

%set up a blocking maze 
blockingMaze = Maze();
blockingMaze.START_STATE = [5,3];
blockingMaze.GOAL_STATE = [0 8];
blockingMaze.oldObstacles = [3 1; 3 2; 3 3; 3 4; 3 5; 3 6; 3 7; 3 8; 3 8];

%new obstacles will block the optimal path 
blockingMaze.newObstacles = [3 2; 3 3; 3 4; 3 5; 3 6; 3 7; 3 8; 3 9;];

%step limit 
blockingMaze.maxSteps =3000;

%when obstacles will change 
blockingMaze.changingPoint = 1000;

%set up parameters for dyna
dynaParams = DynaParams();

dynaParams.alpha = 0.7;
dynaParams.runs = 20; 
dynaParams.timeWeight = 1e-4;%kappa

rewards = zeros(2, blockingMaze.maxSteps);

for run =1:dynaParams.runs
    rewards_ = zeros(2,maxSteps);
    stateActionValues = blockingMaze.stateActionValues; blockingMaze.stateActionValues;
    for i = 1:length(dynaParams.methods)
       blockingMaze.obstacles = blockingMaze.oldObstacles;
       steps =0;
       lastSteps = steps;
       while steps < blockingMaze.maxSteps
          if i ==1
              
          else
              
          end
       end
    end
end






