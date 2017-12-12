/**
 * $Id: GPomdp.cc 127 2007-09-10 17:20:04Z daa $
 */

#include"PGBasics.hh"
#include<time.h>
#include"GPomdp.hh"

namespace libpg {
/**
 * Top level piece of glue for Policy gradient, link the controller
 * and the simulator 
 * @param  Controller 
 * @param  Simulator (problem specific) 
 * @param  discount factor [0, 1]. If 1, trace must be reset
 * occasionally or have 0 mean reward.
 */
GPomdp::GPomdp(Controller* controller, Simulator* simulator, double discount, double stepSize) : RLAlg(controller, simulator, discount, stepSize) {};


/**
 * Do a fixed number of individual steps, aggregating the reward and
 * gradient if learn is true.  No output is done. Could also be termed
 * an "epoch" 
 * @param  returns average rewards for each reward dimension
 * @param  number ofsteps to accumulate gradient for
 *  @param  should learning be turned on, or are we just
 * exploiting/evaluating the current policy?
 */
double GPomdp::doSteps(Vector& totalRewards, int steps, bool learn) {

    Vector action(simulator->getActionDim());
    Vector rewards(simulator->getRewardDim());
    Vector reinforcement(simulator->getRewardDim());

    // Reset for new gradient estimate.
    totalRewards.clear();
    if (learn) controller->resetGrad();
    obs.setSteps(0);

    simulator->getObservation(obs);    
    for (int s = 0; s < steps; s++) {
	obs.setSteps(obs.getSteps() + 1);
	// Must do next line before calling getAction() since
	// getAction will accumulate log action gradients into the
	// trace directly for efficiency.
	if (learn) controller->discountTrace();
	controller->getAction(obs, action, learn);
	simulator->doAction(action); 
	simulator->getReward(rewards);
	totalRewards += rewards;
	// get the new state observation.
	simulator->getObservation(obs);
	if (learn) {
	    reinforcement.assign(rewards - baseline);
	    controller->accumulateGrad(reinforcement, obs);
	}
    }

    totalRewards /= steps;
    return inner_prod(totalRewards, scalar_vector<double>(simulator->getRewardDim(), 1.0))/simulator->getRewardDim();
}



/**
 * Use the GPOMDP estimate to step parameters in the correct
 * direction.
 * @param  steps to estimate gradient for
 * @param  ref to variable counting total learning steps. include
 * any steps used for a line search
 */
double GPomdp::learnCore(int stepsPerEpoch, int& totalSteps) {

    double lastEpochVal;
    Vector totalRewards(simulator->getRewardDim());

    
    // Estimate grad
    lastEpochVal = doSteps(totalRewards, stepsPerEpoch, true);
    totalSteps += stepsPerEpoch;
    
    // Step in right direction
    controller->computeDirection(stepsPerEpoch);
    controller->batchStep();

    if (useAutoBaseline) baseline.assign(totalRewards);

    return lastEpochVal;


}
}
