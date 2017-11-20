classdef TrivialModel < handle
    %trivial model for planning in Dyna-Q
   properties (SetAccess = private)
       model;%a container of a container 
       %for each state, the returned value is another container for the
       %action
       % for each action, we store the next state and reward pairs 
   end
   
   methods 
       %%%%%%%%%%%%%%%%%%%%%Constructor%%%%%%%%%%%%%%%%%%%%%
       function obj = TrivialModel()
           obj.model = containers.Map(); %initialize with empty container 
       end
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
       function feed(obj, currentState, action, newState, reward)
           if isKey(obj.model, num2str(currentState)) == false               
               obj.model(num2str(currentState)) = containers.Map();%initialize this state with empty container               
           end           
           con = obj.model(num2str(currentState));
           con(num2str(action)) = [newState, reward];% this will modify obj.model as well 
       end
       
       function [stateSam, action, newState, reward] =samples(obj)
           stateIndex = randi(numel(keys(obj.model))); %random index for keys in model container
           allKeys = keys(obj.model);
           state = allKeys(stateIndex); 
           stateSam = str2double(state); %convert the string of the state to a numeric 
           state = obj.model(num2str(state));%actually pick a state (which is random)
           
           actionIndex = randi(numel(keys(state)));
           allKeys = keys(state);
           action = allKeys(actionIndex);
           action = state(num2str(action));%actually pick a random action. 
           
           newState = 
           
                      
       end
   end
end