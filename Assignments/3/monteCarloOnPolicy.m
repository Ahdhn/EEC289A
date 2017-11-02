function  [statesUsableAce, statesNoUsableAce] = monteCarloOnPolicy(nEpisodes)
    statesUsableAce = zeros(10,10);
    statesUsableAceCount = ones(10, 10);
    statesNoUsableAce = zeros(10, 10);
    statesNoUsableAceCount = ones(10, 10);
    
    for i=0:nEpisodes
       [state, reward, blah] = play(0,0);
       state(2) = state(2) - 11;
       state(3) = state(3);
       
       if state(1) == 1
           statesUsableAceCount(state(2), state(3)) = statesUsableAceCount(state(2), state(3))+ 1;
           statesUsableAce(state(2), state(3)) = statesUsableAce(state(2), state(3))+ reward;
       else
            statesNoUsableAceCount(state(2), state(3)) = statesNoUsableAceCount(state(2), state(3))+ 1;
            statesNoUsableAce(state(2), state(3)) = statesNoUsableAce(state(2), state(3))+ reward;
       end
    end
    
    statesUsableAce = statesUsableAce ./ statesUsableAceCount;
    statesNoUsableAce = statesNoUsableAce ./statesNoUsableAceCount;
end