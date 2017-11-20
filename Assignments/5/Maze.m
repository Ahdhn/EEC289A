classdef Maze < handle
    properties (SetAccess = public)
        %values 
        WORLD_WIDTH;
        WORLD_HEIGHT;  
        ACTION_UP;  
        ACTION_DOWN;
        ACTION_LEFT;
        ACTION_RIGHT;
        actions;
        START_STATE; 
        GOAL_STATES; 
        obstacles;
        oldObstacles;
        newObstacles;
        changingPoint; %the point where we change the obstacles 
        stateActionValues;%initial state action pair values 
        maxSteps;
        resolution;%track the reslution of the maze
    end
    
    methods 
        %%%%%%%%%%%%%%%%%%%%%Constructor%%%%%%%%%%%%%%%%%%%%% 
        function obj = Maze()
            obj.WORLD_WIDTH =9;
            obj.WORLD_HEIGHT = 6;
            obj.ACTION_UP = 1;
            obj.ACTION_DOWN = 2;
            obj.ACTION_LEFT = 3;
            obj.ACTION_RIGHT = 4;
            obj.actions = [obj.ACTION_UP, obj.ACTION_DOWN,obj.ACTION_LEFT,obj.ACTION_RIGHT];
            obj.START_STATE = [3,1];
            obj.GOAL_STATES = [1,9];
            obj.obstacles = [[2, 3]; [3, 3]; [4, 3]; [1, 8]; [2, 8]; [3, 8]; [5, 6]];                        
            obj.stateActionValues = zeros(obj.WORLD_HEIGHT,obj.WORLD_WIDTH, length(actions));
            obj.maxSteps = inf;
            obj.resolution = 1;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [curState, reward] = takeAction(obj, state, action)
            x = state(1);
            y = state(2);
            if action == obj.ACTION_UP
                x = max(x-1,0);
            elseif action==obj.ACTION_DOWN
                x = min(x + 1, obj.WORLD_HEIGHT - 1);
            elseif action == obj.ACTION_LEFT
                y = max(y - 1, 0);
            elseif action == obj.ACTION_RIGHT
                y = min(y + 1, obj.WORLD_WIDTH - 1);
            end
            
            if ismember([x,y],obj.obstacles)
                x=state(1);
                y=state(2);
            end 
            
            if ismember([x,y],obj.GOAL_STATES)
                reward = 1.0;
            else
                reward = 0.0;
            end            
            curState(1)=x;
            curState(2)=y;
        end
    end
end 