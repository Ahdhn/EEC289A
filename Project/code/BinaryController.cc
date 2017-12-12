#include"PGBasics.hh"
#include"BasicController.hh"
#include"BinaryController.hh"

namespace libpg {

/**
 * @param  parameterised distribution we'll be using today
 */
BinaryController::BinaryController(Approximator* approx) : BasicController(approx) {
    
    // Function approximator must have only 1 output.
    assert(approx->getOutputDim() == 1);
}


/**
 * Implement logistic regression to compute the probability of
 * action=yes, or no.
 * approximator output dimension must = 1
 * @param  observation matrix
 * @param  vector to store action in
 * @param  shall we propogate the log action prob derivative back to
 * approximator?
 */
void BinaryController::getAction(Observation& obs, Vector& action,  bool computeGrad) {
    
    int sampledAct = 0;
    double expDist = 0;
    double pr=0;

    assert(action.size() == 1);
    //Vector tmp(approx->outputs);
    approx->doApprox(obs, dist);

    expDist = exp((dist)[0]);
    pr = 1.0/(expDist + 1.0);
    if (pr >= random()/(double)RAND_MAX) sampledAct = 1;
 
    action[0] = sampledAct;
    
    if (computeGrad) {
	(dist)[0] = pr;
	if (sampledAct) (dist)[0] *= -expDist;
	approx->feedbackGrad(obs, dist);
    }

}


/**
 * For evaluation purposes only. Select the most probable action
 * @param  observation 
 * @param  a vector of length 1 to put in the binary action, 1 or 0.
 */
void BinaryController::getMostProbAction(Observation& obs, Vector& action) { 

    double expDist = 0;
    double pr=0;

    assert(action.size() == 1);
    //Vector tmp(approx->outputs);
    approx->doApprox(obs, dist);

    expDist = exp((dist)[0]);
    pr = 1.0/(expDist + 1.0);
    // If we choose this action with more than prob 0.5, this is the max prob action
    if (pr >= 0.5) action[0] = 1;
    else action[0] = 0;

}
}
