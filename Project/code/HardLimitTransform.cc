/**
 * $Id: RegularizeTransform.cc 51 2007-02-01 21:33:00Z buffet $
 */

#include"PGBasics.hh"
#include"HardLimitTransform.hh"

namespace libpg {

/** 
 * Implements absolute max for parameter values.
 */

HardLimitTransform::HardLimitTransform(Controller* c, double maxp) : TransformController(c) {
    
    this->maxp = maxp;
    this->controller = c;
    params.resize(controller->getNumParams());
    
}


void HardLimitTransform::setMax(double maxp) {
    assert(maxp >= 0.0);
    this->maxp = maxp;
}


void HardLimitTransform::batchStep() {

    // Do the normal step
    controller->batchStep();
 
    checkParams();

}


void HardLimitTransform::instantStep(Vector& rewards) {

    if (norm_1(rewards) > 0.0) {
        controller->instantStep(rewards);
        checkParams();
    }

}


/**
 * The only non-trivial routine.. just make sure all parameter values
 * are within the maximum abosolute val. Limit them if they are.
 */
void HardLimitTransform::checkParams() {
    
    params.clear();

    // Gather up all the parameters
    controller->reduce(params, Approximator::PARAMS);
	    
    bool change = false;

    assert(!UBlasExtras::containsNaN(params));

    // Be careful to maintain the sign
    for (size_t i=0; i < params.size(); i++) {
	if (params[i] > maxp) {
	    params[i] = maxp;
	    change = true;
	}
	else if (params[i] < -maxp) {
	    params[i] = -maxp;
	    change = true;
	}
    }
    
    // Split them out again
    if (change) controller->scatter(params, Approximator::PARAMS);

}
}
