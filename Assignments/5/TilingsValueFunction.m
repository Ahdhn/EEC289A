classdef TilingsValueFunction < handle
    properties (SetAccess = private)
      numOfTilings;
      tileWidth;
      tilingOffset;
      tilingSize;
      params;
      tilings;
   end
   
   methods 
        %%%%%%%%%%%%%%%%%%%%%Constructor%%%%%%%%%%%%%%%%%%%%%
       function obj = TilingsValueFunction(numOfTilings,tileWidth,tilingOffset )
           global N_STATES;
            obj.numOfTilings = numOfTilings;
            obj.tileWidth = tileWidth;
            obj.tilingOffset = tilingOffset;
            obj.tilingSize = N_STATES/ tileWidth+1;
            obj.params = zeros(obj.numOfTilings, obj.tilingSize);
            obj.tilings = -tileWidth+1:tilingOffset:0;           
       end
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
       function stateValue = value(obj, state)
           %get the value of state 
           stateValue =0.0;
           
           %go through all the tilings
           for tilingIndex = 1:length(obj.tilings)
               %find the active tile in current tiling 
               tileIndex = ceil((state -1 - obj.tilings(tilingIndex))/obj.tileWidth);               
               stateValue = stateValue + obj.params(tilingIndex, tileIndex);
           end
       end
       
       function update(obj, delta, state)
           %update parameters
           %each state is covered by some tilings
           %so the delta should be divided equally into each tiling 
           delta = delta / obj.numOfTilings;
           
           %go throught all the tilings 
           for tilingIndex = 1:length(obj.tilings)
               tileIndex = ceil((state -1 - obj.tilings(tilingIndex))/obj.tileWidth);
               obj.params(tilingIndex, tileIndex) = obj.params(tilingIndex, tileIndex) + delta;
           end
       end
   end
end