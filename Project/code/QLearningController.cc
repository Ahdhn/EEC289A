#include"PGBasics.hh"
#include"BasicController.hh"
#include"ValueController.hh"
#include"QLearningController.hh"
#include"eGreedyPolicy.hh"
#include"Sampler.hh"
#include <cfloat>

using namespace std;

namespace libpg {


    /**
     * Q-Learning controller constructor
     * @param approx parameterised distribution we'll be using today
     * @param tdDicount temporal difference discount
     */
    QLearningController::QLearningController(Approximator* approx, Policy* policy, double tdDiscount) : ValueController(approx) {

	this->policy = policy;
	this->tdDiscount = tdDiscount;

	greedyPolicy = new eGreedyPolicy(0.0);
	dist.resize(approx->getOutputDim());
	lastDist.resize(approx->getOutputDim());

    }


    /**
     * Q-Learning controller destructor
     */
    QLearningController::~QLearningController() {
	delete greedyPolicy;
    }


    /**
     * Select the next action based on policy parameter
     * @param obs observation vector
     * @param action single position vector used to return the selected action
     * @param computeGrad flag used to enable learning from experience
     */
    void QLearningController::getAction(Observation& obs, Vector& action,  bool computeGrad) {


	assert(action.size() == 1);
    
	// greedyAct has same dimention as action
	greedyAct.resize(action.size());

	// Learn by experience?
	this->computeGrad = computeGrad;

	// Get state-action values for this state
	approx->doApprox(obs, dist);

	// Get next action according to policy
	policy->getAction(dist, action, obs);

	// Get greedy action
	greedyPolicy->getAction(dist, greedyAct, obs);

	// Temporal difference error calculation
	if (obs.getSteps() <= 1)
	    error = 0.0;
	else {

	    approx->doApprox(lastObs, lastDist);
	    error = lastRew + tdDiscount * dist((int)greedyAct[0]) - lastDist(lastAct);
	}
    
	// Recall last action and observation
	lastAct = (int)action[0];
	lastObs = obs;
    }


    /**
     * Updates the eligibility traces
     */
    void QLearningController::discountTrace() {
	if (computeGrad) {
	    if (lastAct == (int)greedyAct[0]) approx->discountTrace();
	    else approx->resetTrace();
	    dist.clear();
	    dist[lastAct] = 1.0;
	    approx->feedbackGrad(lastObs, dist);
	}
    }


    void QLearningController::instantStep(Vector& reward) {
	if (computeGrad) {
	    lastRew = reward[0];
	    approx->instantStep(error);
	}
    }


    void QLearningController::setDiscount(double lambda) {
	this->lambda = lambda;
	approx->setDiscount(tdDiscount*lambda);
    }


}

