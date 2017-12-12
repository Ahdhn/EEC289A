/**
 * $Id: OLPomdp.cc 127 2007-09-10 17:20:04Z daa $
 */

#include"PGBasics.hh"
#include<time.h>
#include"OLPomdp.hh"

namespace libpg {

/**
 * Top level piece of glue for Policy gradient, link the controller
 * and the simulator 
 * @param  Controller 
 * @param  Simulator (problem specific) 
 * @param  discount factor [0, 1]. If 1, trace must be reset
 * occasionally or have 0 mean reward.
 */
OLPomdp::OLPomdp(Controller* controller, Simulator* simulator, double discount, double stepSize) : RLAlg(controller, simulator, discount, stepSize) {};


/**
 * Do a fixed number of individual steps, aggregating the reward.  No
 * output is done. Could also be termed an "epoch" 
 * @param  number ofsteps to iterate for 
 * @param  should learning be turned on, or are we
 * just exploiting/evaluating the current policy?
 */
double OLPomdp::doSteps(Vector& totalRewards, int steps, bool learn) {


    Vector action(simulator->getActionDim());
    Vector rewards(simulator->getRewardDim());
    Vector reinforcement(simulator->getRewardDim());
    
    totalRewards.clear();


    for (int s = 0; s < steps; s++) {
	// Must do next line before calling getAction() since
	// getAction will accumulate log action gradients into the
	// trace directly for efficiency.
	if (learn && obs.getSteps()>0) controller->discountTrace();
	obs.setSteps(obs.getSteps()+1);
	simulator->getObservation(obs);

	controller->getAction(obs, action, learn);

	if (simulator->doAction(action) != 0) break; 

	// If the rewards are for the curret state, we get the reward
	// for the new state here.  Otherwise, if the reward is state,
	// action, we get the reward for the transition just
	// completed.  Finally, if it's state, action, next state,
	// then we should really move getRewards one line earlier
	// (above doAction).
	// This is a problem I haven't quite figure out how to resolve yet.
	// In algorithms like LSTDQ that assume this latter case,
	// it manually remembers the old reward and stores the
	// reward it's about to get passed.


	simulator->getReward(rewards);
	totalRewards += rewards;
	if (learn) {
	    reinforcement.assign(rewards - baseline);
	    // [daa] Moved the check for 0 reward to the controller level
	    //controller->instantStep(reinforcement, obs);
	    controller->instantStep(rewards);
	}



    }

    totalRewards /= steps;
    return inner_prod(totalRewards, scalar_vector<double>(simulator->getRewardDim(), 1.0))/simulator->getRewardDim();

}

}
