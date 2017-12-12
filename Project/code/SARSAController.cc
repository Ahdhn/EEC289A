#include"PGBasics.hh"
#include"BasicController.hh"
#include"ValueController.hh"
#include"SARSAController.hh"
#include"Sampler.hh"
#include <cfloat>

using namespace std;
 
namespace libpg {


    /**
     * SARSA controller constructor
     * @param approx parameterised distribution we'll be using today
     * @param tdDicount temporal difference discount
     */
    SARSAController::SARSAController(Approximator* approx, Policy* policy, double tdDiscount) : ValueController(approx) {

	this->policy = policy;
	this->tdDiscount = tdDiscount;
	lastRew = 0;
	dist.resize(approx->getOutputDim());
	lastDist.resize(approx->getOutputDim());

    }


    /**
     * SARSA controller destructor
     */
    SARSAController::~SARSAController() {
    }


    /**
     * Select the next action based on policy parameter
     * @param obs observation vector
     * @param action single position vector used to return the selected action
     * @param computeGrad flag used to enable learning from experience
     */
    void SARSAController::getAction(Observation& obs, Vector& action,  bool computeGrad) {

	assert(action.size() == 1);

	// Learn by experience?
	this->computeGrad = computeGrad;


	// Get state-action values for this state
	approx->doApprox(obs, dist);


	// Get next action according to policy
	policy->getAction(dist, action, obs);

	if (computeGrad) {
	    // Temporal difference error calculation
	    if (obs.getSteps() <= 1)
		error = 0.0;
	    else {

		// Get new value estimate of previous state
		approx->doApprox(lastObs, lastDist);
		
		// Compute td error
		error = lastRew + tdDiscount * dist((int)action[0]) - lastDist[lastAct];
		
		//cout<<"lastRew="<<lastRew<<" tddiscount="<<tdDiscount<<" lambda="<<lambda<<" lastAct="<<lastAct<<" action[0]="<<action[0]<<" lastDist="<<lastDist<<" dist="<<dist<<" error="<<error<<endl;
	    }
	}
    
	// Recall last action and observation
	lastAct = (int)action[0];
	lastObs = obs;
    }


    /**
     * Updates the eligibility traces
     */
    void SARSAController::discountTrace() {
	if (computeGrad) {
	    approx->discountTrace();
	    dist.clear();
	    dist[lastAct] = 1.0;
	    approx->feedbackGrad(lastObs, dist);
	}
    }


    void SARSAController::instantStep(Vector& reward) {
	if (computeGrad) {
	    lastRew = reward[0];
	    approx->instantStep(error);
	}
    }

    void SARSAController::setDiscount(double lambda) {
	this->lambda = lambda;
	approx->setDiscount(tdDiscount*lambda);
    }


}
