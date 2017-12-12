/**
 * $Id: RegularizeTransform.cc 127 2007-09-10 17:20:04Z daa $
 */

#include"PGBasics.hh"
#include"RegularizeTransform.hh"

namespace libpg {

/** 
 * Implements Tikhonov style regularisation, either in an online step,
 * or batch step. Regularizer gradient for \theta_i is -q\theta_i
 * Should be able to sit anywhere in the controller stack since it
 * uses the reduce and scatter commands for the params only.
 */

RegularizeTransform::RegularizeTransform(Controller* c, double penalty) : TransformController(c) {
    
    this->penalty = penalty;
    params.resize(controller->getNumParams());
    penaltyVector.resize(controller->getNumParams());
    
}

void RegularizeTransform::setStepSize(double stepSize) {
    // Remember the step size of our children so that
    // we can appropriately scale the regularization.
    this->stepSize = stepSize;
    controller->setStepSize(stepSize);
}

void RegularizeTransform::setPenalty(double penalty) {
    assert(penalty >= 0.0);
    this->penalty = penalty;
}


void RegularizeTransform::batchStep() {

    // Get the pre-step parameters
    if (penalty > 0.0) {
	penaltyVector.clear();
	controller->reduce(penaltyVector, Approximator::PARAMS);
    }

    // Do the normal step
    controller->batchStep();
 
    if (penalty > 0.0) penalise();

}


void RegularizeTransform::instantStep(Vector& rewards) {

    // Get the pre-step parameters
    if (penalty > 0.0) {
	penaltyVector.clear();
	controller->reduce(penaltyVector, Approximator::PARAMS);
    }

    // Do the normal step
    controller->instantStep(rewards);
 
    if (penalty > 0.0) penalise();

}


void RegularizeTransform::penalise() {
    
    // Now post-process the new parameters
    // Scale to alpha and penalty
    penaltyVector *= -stepSize*penalty;

    // Do the penalisation.
    params.clear();
    controller->reduce(params, Approximator::PARAMS);
    params += penaltyVector;
    controller->scatter(params, Approximator::PARAMS);

}
}
