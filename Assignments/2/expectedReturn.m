%state = [number of cars in 1st location, number of cars in 2nd location]
%action = +ve if moving cars from 1st to 2nd locations 
%         -ve if moving cars from 2nd to 1st locations 
%stateValue is state value array 

function returns = expectedReturn(state,action,stateValue)
                                 
    global poissonBackup MOVE_CAR_COST POISSON_UP_BOUND RENTAL_CREDIT...
           RENTAL_REQUEST_FIRST_LOC RENTAL_REQUEST_SECOND_LOC DISCOUNT...
           MAX_CARS RETURNS_FIRST_LOC RETURNS_SECOND_LOC SECOND_PARKING_COST;
    
    returns =0.0;%initialize
    returns = returns - MOVE_CAR_COST *abs(action);%cost of moving cars
    %since there is someone can move the car for free, then we shouln't
    %discount it, only if we are actually moving any cars from 1st to 2nd
    %location (which is given positive action value)
    if action>0
       returns = returns + MOVE_CAR_COST;
    end
    
    %all possible rental requests
    for rentalRequestFirstLoc = 0:POISSON_UP_BOUND-1
        for rentalRequestSecondLoc  = 0:POISSON_UP_BOUND-1
            %moving cars
            numOfCarsFirstLoc  = min(state(1) - action,MAX_CARS);
            numOfCarsSecondLoc = min(state(2) + action, MAX_CARS);
            
            %rental requests should be less than actual number of cars
            realRentalFirstLoc = min(numOfCarsFirstLoc, rentalRequestFirstLoc);
            realRentalSecondLoc = min(numOfCarsSecondLoc, rentalRequestSecondLoc);
                        
            numOfCarsFirstLoc =  numOfCarsFirstLoc - realRentalFirstLoc;
            numOfCarsSecondLoc =  numOfCarsSecondLoc - realRentalSecondLoc;
            
            %credits for renting 
            reward = (realRentalFirstLoc + realRentalSecondLoc)*RENTAL_CREDIT;
            
            %punish the configuration that keep more than 10 cars in any of
            %the two locations 
            
            if numOfCarsFirstLoc > 10
                reward = reward -SECOND_PARKING_COST;
            end
            if numOfCarsSecondLoc > 10
                reward = reward -SECOND_PARKING_COST;
            end
            
            %probability for current combination of rental request
            prob = poisson(rentalRequestFirstLoc, RENTAL_REQUEST_FIRST_LOC)*...
                   poisson(rentalRequestSecondLoc, RENTAL_REQUEST_SECOND_LOC);
               %get returned cars, those cars can be used for renting tomorrow
               returnedCarsFirstLoc= RETURNS_FIRST_LOC;
               returnedCarsSecondLoc = RETURNS_SECOND_LOC;
               numOfCarsFirstLoc = 1 + min(numOfCarsFirstLoc + returnedCarsFirstLoc, MAX_CARS);
               numOfCarsSecondLoc = 1 + min(numOfCarsSecondLoc + returnedCarsSecondLoc, MAX_CARS);
               returns = returns + prob*(reward + DISCOUNT*stateValue(numOfCarsFirstLoc, numOfCarsSecondLoc));                                    
        end        
    end
end