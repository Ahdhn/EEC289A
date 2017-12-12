
// $Id: TransformController.cc 127 2007-09-10 17:20:04Z daa $


#include"PGBasics.hh"

using namespace std;

namespace libpg {

    /**
     * @param  underlying controller.
     */
    TransformController::TransformController(Controller* controller) {
	this->controller = controller;
    }


    /**
     * Just pass to underlying controller
     */
    void TransformController::getAction(Observation& obs, 
					Vector& action,  
					bool computeGrad) {
	controller->getAction(obs, action, computeGrad);
    }

    void TransformController::computePartialsAndFeedback(Observation& obs, int sampledAct, double weight) {
	controller->computePartialsAndFeedback(obs, sampledAct, weight);
    }

    /**
     * Just pass to underlying controller
     * @param  Observation vector
     * @param  empty vector to load with the most probably action
     */
    void TransformController::getMostProbAction(Observation& obs, Vector& action) { 
	controller->getMostProbAction(obs, action);
    }


    /**
     * Just pass to underlying controller
     */
    int TransformController::getInputDim() {
	return controller->getInputDim();
    }

    /**
     *
     * The remaining methods just pass onto equiv Controller
     * functions. We have to do it in this ugly way because more advanced
     * multiple agent controllers will need to override these functions to
     * update all agents in turn.
     *
     */

    void TransformController::discountTrace() {
	controller->discountTrace();
    }


    void TransformController::setDiscount(double discount) {
	controller->setDiscount(discount);
    }


    void TransformController::accumulateGrad(Vector& rewards, Observation& newObs) {
	controller->accumulateGrad(rewards, newObs);
    }


    void TransformController::resetGrad() {
	controller->resetGrad();
    }


    void TransformController::instantStep(Vector& rewards) {
	controller->instantStep(rewards);
    }


    void TransformController::computeDirection(int steps) {
	controller->computeDirection(steps);
    }


    void TransformController::batchStep() { 
	controller->batchStep(); 
    }


    void TransformController::setStepSize(double stepSize) { 
	controller->setStepSize(stepSize); 
    }


    void TransformController::resetTrace() { 
	controller->resetTrace(); 
    }


    void TransformController::resetParams() { 
	controller->resetParams(); 
    }


    void TransformController::randomizeParams(double maxRand) { 
	controller->randomizeParams(maxRand); 
    }


    double TransformController::getMaxParam() {
	return controller->getMaxParam();
    }


    int TransformController::getNumParams() {
	return controller->getNumParams();
    }


    void TransformController::write(ostream& o) {
	controller->write(o);
    }


    void TransformController::read(istream& o) {
	controller->read(o);
    }


    void TransformController::reduce(Vector& v, Approximator::StatsEnum s) {
	controller->reduce(v, s);
    }


    void TransformController::scatter(Vector& v, Approximator::StatsEnum s) {
	controller->scatter(v, s);
    }

}
