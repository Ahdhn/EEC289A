classdef DynaParams < handle
    %all dyna parameters
    properties (SetAccess = public)
        gamma;%discount
        epsilon; %propability for exploration
        alpha;%step size
        timwWeight;%wieght of elapsed time 
        planningSteps; %n-step planning
        runs;%avergae over several independent runs 
        method; %algorithm name 
        theta;%threshold for priority queue        
    end
    
    methods 
        %%%%%%%%%%%%%%%%%%%%%Constructor%%%%%%%%%%%%%%%%%%%%% 
        function obj = DynaParams()
            obj.gamma = 0.95;
            obj.epsilon = 0.1;
            obj.alpha = 0.1;
            obj.timwWeight = 0; 
            obj.planningSteps = 5;
            obj.runs = 10; 
            obj.method = ['Dyna-Q','Dyna-Q+'] ; 
            obj.theta = 0;   
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    
end