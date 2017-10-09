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
    end

    methods
       %Constructor%
       function obj = Bandit(kArm,stepSize,sampleAverages,eps, initial)
           obj.kArm =kArm;
           obj.stepSize=stepSize;
           obj.indices = 1:kArm;                         
           obj.sampleAverages = sampleAverages;
           obj.time =0;
           obj.averageReward =0;
           obj.trueReward =0;
           obj.qTrue =[];
           obj.qEst = initial*ones(1,kArm);
           obj.eps=eps;
           obj.actionCount = zeros(1,kArm);
           
           for i=1:kArm
               obj.qTrue(i) = normrnd(0,1) + obj.trueReward;
           end
           
           obj.bestAction = max(obj.qTrue);
       end
       
       function getAction(obj)
           %explore 
           if obj.eps >0
               
           end
          
       end
       
    end
    
    
end