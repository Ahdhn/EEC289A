#include"PGBasics.hh"
#include"Sampler.hh"
#include"RegressionController.hh"


using namespace std;
namespace libpg {

    /**
     * @param  approximators vector of Approximator* 
     * @param  if true, expect as many colums in the observations as there are approximators
     * and each controller only gets its col as an observation.
     * If false, all approximators get the same observation.
     */
    RegressionController::RegressionController(Approximator* approx) {

	this->approx = approx;


    }


    /**
     * Here we just compute the function approximator output
     */
    void RegressionController::getAction(Observation& obs, Vector& action,  bool computeGrad) {
    
	approx->doApprox(obs, action);
	lastObs = obs;
    
    }
  
    /**
     * This is a bit of a hack to allow us to call
     * computePartialsAndFeedback() with a regression controller.
     * Since it's not in the spec for a Controller, you have to cast
     * Controllers to RegressionControllers in order to use it. I
     * should think of a neater way to do this (probably changing int
     * sampledAct to a vector.
     */
    void RegressionController::setLastReward(Vector& lastRewards) {
	this->lastRewards = lastRewards;
    }

    /**
     * This is a bit odd in the context of regression. The action
     * isn't really how we compute partials, it's the reward (which is
     * treated directly as vector of partials). So what it does is computes
     * the partial on the basis of the last reward that was specified via 
     */
    void RegressionController::computePartialsAndFeedback(Observation& obs, int sampledAct, double weight) {

	resetTrace();
	approx->feedbackGrad(obs, lastRewards);
    
    }
    
    
    /**
     * Pick the controller that has the highest probability of being chosen.
     * @param  observation vector
     * @param  empty vector length 1 to put the action in.
     */
    void RegressionController::getMostProbAction(Observation& obs, Vector& action) {

	getAction(obs, action, false);
    
    }
    
    
    int RegressionController::getInputDim() {
	return approx->getInputDim();
    }
   

    /**
     *
     * The remaining methods just pass onto equiv Approximator
     * functions. We have to do it in this ugly way because more advanced
     * multiple agent approximators will need to override these functions to
     * update all agents in turn.
     *
     */

    void RegressionController::discountTrace() {
    }


    void RegressionController::setDiscount(double discount) {
    }


    void RegressionController::accumulateGrad(Vector& rewards, Observation& newObs) {

	resetTrace();
	approx->feedbackGrad(lastObs, rewards);
	approx->accumulateGrad(1.0, newObs);	
    }


    void RegressionController::resetGrad() {
	approx->resetGrad();
    }


    void RegressionController::instantStep(Vector& rewards) {
	resetTrace();
	approx->feedbackGrad(lastObs, rewards);
	approx->instantStep(1.0);
    }


    void RegressionController::computeDirection(int steps) { 
	approx->computeDirection(steps);
    }


    void RegressionController::batchStep() { 
	approx->batchStep();
    }

    void RegressionController::setStepSize(double stepSize) { 
	approx->setStepSize(stepSize);
    }


    void RegressionController::resetTrace() { 
	approx->resetTrace();
    }


    void RegressionController::resetParams() { 
	approx->resetParams();
    }
    

    void RegressionController::randomizeParams(double maxRand) { 
	approx->randomizeParams(maxRand);
    }
    

    double RegressionController::getMaxParam() {
	return approx->getMaxParam();
    }


    int RegressionController::getNumParams() { 	
	return approx->getNumParams();
    }


    void RegressionController::write(ostream& o) {
	approx->write(o);
    }
  

    void RegressionController::read(istream& o) {
	approx->read(o);
    }  

    /**
     * @param s what to set, can't do GRADS here.
     */
    void RegressionController::scatter(Vector& v, Approximator::StatsEnum s) {
	approx->scatter(v, s);
    }

    void RegressionController::reduce(Vector& v, Approximator::StatsEnum s) {
	approx->reduce(v, s);
    }

}
