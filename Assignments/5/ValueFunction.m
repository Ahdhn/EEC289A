classdef ValueFunction < handle
    properties (SetAccess = private)
       numOfGroups;
       groupSize;
       params;
   end
   
   methods 
        %%%%%%%%%%%%%%%%%%%%%Constructor%%%%%%%%%%%%%%%%%%%%%
       function obj = ValueFunction(numOfGroups)
           global N_STATES
           obj.numOfGroups = numOfGroups;
           obj.groupSize = N_STATES / numOfGroups;
           obj.params = zeros(1,numOfGroups); %theta     
       end
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function val = value(obj, state)
           %get the value of state 
           global END_STATES
           if ismember(state, END_STATES)
               val =0;
           else               
               groupIndex = ceil((state-1) / obj.groupSize);                
               val = obj.params(groupIndex);
           end           
       end
       
       function update(obj, delta, state)
           %update parameters            
           groupIndex = ceil((state -1)/ obj.groupSize);
           
           obj.params(groupIndex) = obj.params(groupIndex) + delta;
       end
   end
end