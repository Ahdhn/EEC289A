function [state, reward, playerTrajectory] = play(useBehaviourPolicy, initialState )
    global policyDealer;
    playerSum = 0; %sum of player
    usableAcePlayer = false;%does the player uses Ace as 11    
    playerTrajectory =inf; %to indicate we have not appended anything here
    
    
    %dealer status  
    dealerCard1 =0;
    dealerCard2 =0;
    usableAceDealer = false;
    
    getCard=@()(min(randi(14),10));
    
    if initialState == 0
        %generate random initial state
        numOfAce =0;
        
        %initialize cards of player
        while playerSum < 12
            %if sum of player is less than 12, always hit
            card = getCard();
            
            if card == 1%if an ace, use it as 11
                numOfAce = numOfAce+1;
                card =11;
                usableAcePlayer = true;
            end
            playerSum = playerSum+card;
        end
        
        if playerSum > 21
            %if player's sum is larger than 21, he must hold at least one
            %Ace, two Aces are possible
            playerSum = playerSum -10;
            
            if numOfAce == 1%if the player has only one Ace, then s/he 
                %have usable Ace anymore
                usableAcePlayer = false;
            end
        end
        
        %initialize cards for dealer, suppose dealer will show the first
        %card he gets
        dealerCard1 = getCard();
        dealerCard2 = getCard();
        
    else
        usableAcePlayer = initialState(1);
        playerSum = initialState(2);
        dealerCard1 = initialState(3);
        dealerCard2 = getCard();
          
    end
    
    %initial state of the game 
    state = [usableAcePlayer, playerSum, dealerCard1];
    
    
    %initialize dealer's sum
    dealerSum = 0;
    if dealerCard1 == 1 && dealerCard2 ~= 1
        dealerSum = dealerSum + 11 + dealerCard2;
        usableAceDealer = true;
    elseif dealerCard1 ~= 1 && dealerCard2 == 1
        dealerSum = dealerSum + dealerCard1 + 11;
        usableAceDealer = true;
    elseif dealerCard1 == 1 && dealerCard2 == 1
        dealerSum = dealerSum + 1 + 11;
        usableAceDealer = true;
    else
        dealerSum = dealerSum + dealerCard1 + dealerCard2;
    end
      
    %%%% GAME STARTS
    
    %player's turn
    while true                
        if useBehaviourPolicy == 0
            %get action based on current sum
            action = targetPolicyPlayer(playerSum);
        else
            action = behaviorPolicyPlayer();
        end
        
        %update trajectory for importance sampling
        if playerTrajectory == inf
            playerTrajectory = [action,usableAcePlayer,playerSum,dealerCard1];
        else
            playerTrajectory = [playerTrajectory; [action,usableAcePlayer,playerSum,dealerCard1]];
        end
        
        if action == 1
            break;
        end
        
       playerSum = playerSum + getCard(); 
       
       %player busts
       if playerSum >21
           %if player has a usable Ace, use it as 1 to avoid busting and
           %continue
           if usableAcePlayer == true
              playerSum = playerSum -10;
              usableAcePlayer = false;
           else
               %otherwise player loses
               reward = -1;
               return;
           end
       end       
    end
    
    %dealer turn 
    while true
        %get action based on current sum
        action = policyDealer(dealerSum);
        if action == 1
            break;
        end
        
        %if hit, get a new card
        dealerSum = dealerSum + getCard();
        
        %dealer busts
        if dealerSum > 21
            if usableAceDealer == true
                %if dealer has a usable Ace, use it as 1 to avoid busting
                %and continue
                dealerSum = dealerSum -10;
                usableAceDealer =false;
            else
                %otherwise dealer loss
                reward =1;                
                return;
            end
        end
    end 
    
    %compare the sum between the player and the dealer
    if playerSum > dealerSum
        reward = 1;
    elseif playerSum == dealerSum
        reward =0;
    else
        reward =-1;
    end
    
end

function target = targetPolicyPlayer(playerSum)
    global policyPlayer;
    target = policyPlayer(playerSum);
end
function action = behaviorPolicyPlayer()
    if nbinrnd(1,0.5) == 1
        action = 1;
    else
        action =0;
    end 
    
end