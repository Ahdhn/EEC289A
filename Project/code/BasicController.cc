#include"PGBasics.hh"
#include"BasicController.hh"
#include"Sampler.hh"


using namespace std;
namespace libpg {

    /**
     * @param  parameterised distribution we'll be using today
     */
    BasicController::BasicController(Approximator* approx) {

	this->approx = approx;
	dist.resize(approx->getOutputDim());
	distribution.resize(approx->getOutputDim());
    
    }


    /**
     * Where the real coding is done to implement softmax
     */
    void BasicController::getAction(Observation& obs, Vector& action,  bool computeGrad) {
    
	int sampledAct;
	double oneNorm = 0;

	assert(action.size() == 1);
	approx->doApprox(obs, dist);

 
	// Element wise exponentiation and 1-norm, first step of softmax
	int i, imax=dist.size();
	for (i = 0; i < imax; i++) {
	    if (obs.getEligible()[i] != 0) {
		dist[i] = exp(dist[i]);
		oneNorm += dist[i];
	    } else {
		dist[i] = 0.0;
	    }
	}


	// Pick an action. This loop checks for the rare circumstance that an non-eligible action is
	// chosen in a FactoredApproximator
	do {
	    sampledAct  = Sampler::discrete(dist, oneNorm);
	} while (obs.getEligible()[sampledAct] == 0.0);

	assert(sampledAct < approx->getOutputDim());
	action[0] = sampledAct;

	dist /= oneNorm;

	/*cerr<<"BasicController::dist="<<dist<<endl;*/
	distribution = dist;

	if (computeGrad) computePartialsAndFeedback(obs, sampledAct);

    }


    void BasicController::computePartialsAndFeedback(Observation& obs, int sampledAct, double weight) {

	dist *= -1.0;
	//noalias(dist) = -dist;
	dist[sampledAct] += 1.0;
	// Gradient will be added directly to trace
	approx->feedbackGrad(obs, dist);

    }


    /**
     * Return the most probable action for evaluation purposes only.
     * @param  Observation vector
     * @param  empty vector to load with the most probably action
     */
    void BasicController::getMostProbAction(Observation& obs, Vector& action) { 

	assert(action.size() == 1);
	approx->doApprox(obs, dist);

	// Pick most probable action.
	action[0] = 0;

	for (size_t i = 1; i < dist.size(); i++) {
	    if (dist[i] > dist[(int)action[0]]) action[0] = i;
	}
       
    }


    /**
     * Just return the input dimensionality of the approximator
     */
    int BasicController::getInputDim() {
	return approx->getInputDim();
    }


    /**
     * The remaining methods just pass onto equiv Approximator
     * functions. We have to do it in this ugly way because more advanced
     * multiple agent controllers will need to override these functions to
     * update all agents in turn.
     *
     */

    void BasicController::discountTrace() {
	approx->discountTrace();
    }

    void BasicController::setDiscount(double discount) {
	approx->setDiscount(discount);
    }

    void BasicController::accumulateGrad(Vector& rewards, Observation& newObs) {
	assert(rewards.size() == 1);
	approx->accumulateGrad(rewards[0], newObs);
    }

    void BasicController::resetGrad() {
	approx->resetGrad();
    }

    void BasicController::instantStep(Vector& rewards) {
	assert(rewards.size() == 1);
    
	if (rewards[0] != 0.0) approx->instantStep(rewards[0]);
    }


    void BasicController::computeDirection(int steps) {
	approx->computeDirection(steps);
    }

    void BasicController::batchStep() { 
	approx->batchStep(); 
    }

    void BasicController::setStepSize(double stepSize) { 
	approx->setStepSize(stepSize); 
    }


    void BasicController::resetTrace() { 
	approx->resetTrace(); 
    }

    void BasicController::resetParams() { 
	approx->resetParams(); 
    }


    void BasicController::randomizeParams(double maxRand) { 
	approx->randomizeParams(maxRand); 
    }


    double BasicController::getMaxParam() {

	return approx->getMaxParam();

    }

    int BasicController::getNumParams() {
	return approx->getNumParams();
    }


    void BasicController::write(ostream& o) {
	approx->write(o);
    }

    void BasicController::read(istream& o) {
	approx->read(o);
    }

    void BasicController::reduce(Vector& v, Approximator::StatsEnum s) {

	switch (s) {
	case Approximator::DIST:
            v = distribution;
            break;
        default:
	    approx->reduce(v, s);
	}

    }

    void BasicController::scatter(Vector& v, Approximator::StatsEnum s) {
	approx->scatter(v, s);
    }


}
