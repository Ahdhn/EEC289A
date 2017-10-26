%%%%%%%ECE 289A - An Introduction to Reinforcement Learning
%%%HW 2 by Ahmed H. Mahmoud 

%%Re-generating Figure 3.8 by converting a simplified version of the code in 
%https://github.com/ShangtongZhang/reinforcement-learning-an-introduction/blob/master/chapter03/GridWorld.py
%into matlab code 

clc;
clear;
close all;

GRID_SIZE = 5;
A_POS = [1 2];
A_PRIME_POS = [5 2];
B_POS = [1 4];
B_PRIME_POS = [3 4];
discount = 0.9;
actions = ["L" "U" "R" "D"];

%initialize the grid 
for x = 1:GRID_SIZE  
    for y = 1:GRID_SIZE
        actionProb(x,y).L = 0.25;
        actionProb(x,y).U = 0.25;
        actionProb(x,y).R = 0.25;
        actionProb(x,y).D = 0.25;
    end
end

%initialize the rewards and the next states
for x = 1:GRID_SIZE
    for y = 1:GRID_SIZE        
       
        if x == 1
            nextState(x,y).U = [x,y];
            actionReward(x,y).U = -1.0;            
        else
            nextState(x,y).U = [x-1,y];
            actionReward(x,y).U = 0.0; 
        end

        if x == GRID_SIZE
             nextState(x,y).D = [x,y];
            actionReward(x,y).D = -1.0;             
        else
             nextState(x,y).D = [x+1,y];
             actionReward(x,y).D = 0.0;             
        end

        if y == 1
            nextState(x,y).L = [x,y];
            actionReward(x,y).L = -1.0;             
        else
            nextState(x,y).L = [x,y-1];
            actionReward(x,y).L = 0.0;             
        end

        if y == GRID_SIZE
            nextState(x,y).R = [x,y];
            actionReward(x,y).R = -1.0;             
        else
            nextState(x,y).R = [x,y+1];
            actionReward(x,y).R = 0.0;             
        end

        if A_POS(1)==x && A_POS(2) == y
            nextState(x,y).L = A_PRIME_POS;
            nextState(x,y).R = A_PRIME_POS;
            nextState(x,y).D = A_PRIME_POS;
            nextState(x,y).U = A_PRIME_POS;
            actionReward(x,y).L = 10.0;  
            actionReward(x,y).R = 10.0;
            actionReward(x,y).D = 10.0;
            actionReward(x,y).U = 10.0;            
        end

        if B_POS(1)==x && B_POS(2) == y
            nextState(x,y).L = B_PRIME_POS;
            nextState(x,y).R = B_PRIME_POS;
            nextState(x,y).D = B_PRIME_POS;
            nextState(x,y).U = A_PRIME_POS;
            actionReward(x,y).L = 5.0;  
            actionReward(x,y).R = 5.0;
            actionReward(x,y).D = 5.0;
            actionReward(x,y).U = 5.0;            
        end
        
    end
end


%game loop
myGrid = zeros(GRID_SIZE);
while true
   newGrid = zeros(GRID_SIZE);
   
   for x=1:GRID_SIZE
       for y=1:GRID_SIZE
           value=zeros(1,4);
           
           %for the four possible actions
           newPos = nextState(x,y).L;
           value(1)=actionReward(x,y).L + discount*myGrid(newPos(1),newPos(2));
           
           newPos = nextState(x,y).R;
           value(2)=actionReward(x,y).R + discount*myGrid(newPos(1),newPos(2));
           
           newPos = nextState(x,y).D;
           value(3)=actionReward(x,y).D + discount*myGrid(newPos(1),newPos(2));
           
           newPos = nextState(x,y).U;
           value(4)=actionReward(x,y).U + discount*myGrid(newPos(1),newPos(2));
           
           newGrid(x,y)= max(value);          
       end
   end
   
   %check for convergence 
   if isequal(myGrid,newGrid)
       disp('Optimal Solutions to the gridworld example:');
       disp(myGrid);
       break;
   end   
   %copy the new state 
   myGrid = newGrid;   
end