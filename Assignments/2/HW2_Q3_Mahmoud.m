%%%%%%%ECE 289A - An Introduction to Reinforcement Learning
%%%HW 2 by Ahmed H. Mahmoud 

%%Solving Exercise 4.4 by following the steps in 
%https://github.com/ShangtongZhang/reinforcement-learning-an-introduction/blob/master/chapter04/CarRental.py
%and adding the new conditions

clc;
clear;
close all;

global MAX_CARS RENTAL_REQUEST_FIRST_LOC RENTAL_REQUEST_SECOND_LOC...
    RETURNS_FIRST_LOC RETURNS_SECOND_LOC DISCOUNT RENTAL_CREDIT...
    MOVE_CAR_COST poissonBackup POISSON_UP_BOUND SECOND_PARKING_COST;

MAX_CARS = 20;%maximum number of cars in each location
MAX_MOVE_OF_CARS = 5;%maximum number of cars to move during night
RENTAL_REQUEST_FIRST_LOC = 3;%expectation for rental requests in 1st location
RENTAL_REQUEST_SECOND_LOC = 4;%expectation for rental requests in 2nd location
RETURNS_FIRST_LOC = 3;%expectation for number of cars returned in 1st location
RETURNS_SECOND_LOC = 2;%expectation for number of cars returned in second location
DISCOUNT = 0.9;%discount rate
RENTAL_CREDIT = 10;%credit earned by a car
MOVE_CAR_COST = 2;%cost of moving a car
SECOND_PARKING_COST =4;%the cost fo using the second parking space 

policy = zeros(MAX_CARS + 1,MAX_CARS + 1);%current policy 
stateValue = zeros(MAX_CARS+1,MAX_CARS+1);%current state value
actions = -MAX_MOVE_OF_CARS:MAX_MOVE_OF_CARS;%all possible actions 

POISSON_UP_BOUND = 11;%bound on poisson dist; if n is greater than this, then the probability if getting n is truncated to 0

poissonBackup = containers.Map();%Probability for poisson distribution

%all possible states 
numStates =0;
for x = 0:MAX_CARS
    for y = 0:MAX_CARS
        numStates = numStates+1;
        states(numStates,1) = x;
        states(numStates,2) = y;
    end
end

newStateValue = zeros(MAX_CARS+1,MAX_CARS+1);
improvePolicy = false;
policyImprovementInd = 0;
numEval=0;
pi_num=0;
plotPolicy(policy, pi_num);
 
while true
    if improvePolicy == true
        %improve the policy       
        disp(['Policy improvement ', num2str(policyImprovementInd)]);
        policyImprovementInd = policyImprovementInd +1;
        newPolicy = zeros(MAX_CARS+1,MAX_CARS+1);
        
        for stateID = 1:numStates
            x = states(stateID,1);
            y = states(stateID,2);
            
            %go through all actions and select the best one 
            myBestActionReturn = -inf;
            myBestAction = -inf; 
            for ac = 1: length(actions)
                action = actions(ac);
                if (action >=0 && x>=action ) || (action <0 && y >=abs(action))
                    thisActionReturn= expectedReturn([x,y],action,stateValue);   
                    if thisActionReturn > myBestActionReturn
                        myBestActionReturn = thisActionReturn;
                        myBestAction = action;
                    end
                end                
            end            
            newPolicy(x+1,y+1) = myBestAction;            
        end
        %if policy is stable
        policyChanges=sum(sum(newPolicy ~=policy)); 
        if policyChanges == 0
            policy = newPolicy;
            pi_num= pi_num+1;      
            plotPolicy(policy, pi_num);
            break;
        end
        
        policy = newPolicy;        
        disp(['Policy for ',num2str(policyChanges), ' states changed']);
                      
        pi_num= pi_num+1;         
        plotPolicy(policy, pi_num);    
        
        improvePolicy = false;
    end  
    
    
    %start policy evaluation
    for stateID = 1:numStates
        x = states(stateID,1);
        y = states(stateID,2);
        newStateValue(x+1,y+1) = expectedReturn([x y],policy(x+1,y+1), stateValue);   
    end
    if sum(abs(newStateValue - stateValue)) < 1e-4
        stateValue = newStateValue;
        improvePolicy = true;
        continue;
    end      
    numEval = numEval+1;
    
    stateValue = newStateValue;      
end
contour3(stateValue,1000);