classdef Bandit < handle
    
    properties (SetAccess = private)
        %values
       kArm; %number of arms 
       stepSize; %time step size 
       indices;%arm id
       sampleAverages; %true if using sample averages (not constant step)
       time;
       averageReward; 
       trueReward; 
       qTrue;    %real reward for each action 
       qEst;     %estimation of each action
       eps; %probablity of exploration 
       initial; %initial estimation for each action 
       actionCount; %number of times an action has been chosen
       bestAction;
       UCB_param;%if <0, then don't use UCB
       isGradient;%use the graident method 
    end

    methods
       %%%%%%%%%%%%%%%%%%%%%Constructor%%%%%%%%%%%%%%%%%%%%%
       function obj = Bandit(kArm,stepSize,sampleAverages,eps, initial,UCB_param, isGradient)
           obj.kArm =kArm;
           obj.stepSize=stepSize;
           obj.indices = 1:kArm;                         
           obj.sampleAverages = sampleAverages;
           obj.time =0;
           obj.averageReward =0;
           obj.trueReward =0;
           obj.qTrue =[];
           obj.qEst = initial*ones(1,kArm);
           obj.eps=eps*eps;
           obj.actionCount = zeros(1,kArm);
           obj.UCB_param = UCB_param;
           for i=1:kArm
               obj.qTrue(i) = normrnd(0,1) + obj.trueReward;
           end
           
           obj.bestAction = max(obj.qTrue);
           obj.isGradient = isGradient;
       end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
        
        
       %%%%% Which action to take (which arm to pull)%%%%%%%%
       function id = getAction(obj)
           %explore 
           id = -1;
           if obj.eps >0 % if we are working with eps-greedy
               %check if we should explore
               if nbinrnd(1,obj.eps) ==1
                   rand_shuffle = randperm(length(obj.indices));
                   id = rand_shuffle(1);                   
               end                   
           end
           
           %exploit
           if id <0 && obj.UCB_param >=0 %is it UCB 
               UCB_est = obj.qEst + obj.UCB_param * sqrt( log(obj.time+1)./(obj.actionCount+1));
               [maxval,id] = max(UCB_est);
           end  
           
           % TODO add the gradient method 
           
           if id <0 %if it is not UCB and we did not explore
               [maxval,id] = max(obj.qEst);               
           end
       end
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
       
       
       
       %%take an action, update estimation for this action%%%
       function reward= takeAction(obj, action)
           reward = randn() + obj.qTrue(action); %rand reward from normal distribution
           obj.time = obj.time+1; %increment the time 
           obj.averageReward = ((obj.time-1.0)/obj.time)*obj.averageReward + reward/obj.time; %compute the current avergare reward
           obj.actionCount(action) =obj.actionCount(action) +1; %increment this action count
           
           if obj.sampleAverages == true
               %update the estimation using sample averages 
               obj.qEst(action) = obj.qEst(action) + (  (1.0/obj.actionCount(action)) * (reward-obj.qEst(action)));   
               
               %TODO add the gradient method here
               
           else
               %update the estimation with constant step size
               obj.qEst(action) = obj.qEst(action) + (obj.stepSize*(reward-obj.qEst(action))); 
           end          
       end       
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
    end
    
    
end