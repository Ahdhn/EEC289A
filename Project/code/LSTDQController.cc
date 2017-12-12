/**
 * $Id: LSTDQController.cc 68 2007-02-12 06:26:33Z daa $
 */

#include"PGBasics.hh"
#include"LSTDQController.hh"

// An initial value for the AInverse matrix corresponding to an initial A
// which is a small fraction of the identity matrix 
#define IDENT_SCALE (1.0/0.001)

using namespace std;
namespace libpg {

    /**
     * Constructor. Set up additional vectors and matrices needed to
     * implement natural actor critic.  
     * @param actions the
     * LSTDQCOntroller needs to know how many actions it has since the
     * approximator doesn't know about actions in this case, only
     * Q-values
     */
    LSTDQController::LSTDQController(Approximator* approx, Policy* policy, double tdDiscount) : ValueController(approx) {

	this->policy = policy;
	this->tdDiscount = tdDiscount;

	// The total size of the param vector should be actions*getNumInputs()
	// in the linear case
	lstdqDimension = approx->getNumParams();
	
	I.resize(lstdqDimension, lstdqDimension);
	I.assign(identity_matrix<double>(lstdqDimension));
	I *= IDENT_SCALE; 

	//A.resize(lstdqDimension, lstdqDimension);
	AInverse.resize(lstdqDimension, lstdqDimension);
	AInverse.clear();

	values.resize(approx->getOutputDim()); // The number of actions
	
	lstdqTD.resize(lstdqDimension);
	lstdqTrace.resize(lstdqDimension);
	lstdqTrace.clear();
	bVec.resize(lstdqDimension);
	bVec.clear();
	lstdqWeights.resize(lstdqDimension);
	lstdqWeights.clear();

	previousPhi.resize(lstdqDimension);
	previousPhi.clear();
	currentPhi.resize(lstdqDimension);
	currentPhi.clear();
	
	instantStepRequired = false;
	lstdqSteps = 0;
		
    }
    
    
    /**
     * LSTDQ delays it's discounting till the end of a updateStatistics()
     */
    void LSTDQController::discountTrace() {
	
    }

    /**
     * Just a wrapper that intercepts the new observation and completes
     * the instant step from the last round if necessary.
     */
    void LSTDQController::getAction(Observation& obs, Vector& action, bool learn) {

	assert(obs.getFeatures().size1() == (size_t)getInputDim());
	
	// Doing an approx step computes values
	approx->doApprox(obs, values);


	// Now values stores the q values over all the actions for the
	// current state (features)
	policy->getAction(values, action, obs);


	if (learn) {
	
	    // Compute the basis feature phi by a
	    // feedbackGrad(). Relies on this being a linear
	    // approximator!!
	    values.clear();
	    values[(size_t)action[0]] = 1.0;

	    // feedbackGrad puts the result in the trace
	    approx->resetTrace();
	    approx->feedbackGrad(obs, values);
	    
	    // Save the previous phi 
	    previousPhi = currentPhi;

	    // extract the phi for the basis feature and action just chosen.
	    currentPhi.clear();
	    approx->reduce(currentPhi, Approximator::TRACES);

	    // Do a step if required
	    if (instantStepRequired) delayedInstantStep(lastReward, obs.getSteps());
	    
	}
	
    }
    

    
    // Reset just the eligibility trace. Note that LSTDQ typically
    // assumes \lambda=0, so this would do nothing in that case.
    void LSTDQController::resetTrace() {
	lstdqTrace.clear();
	AInverse.assign(I); 
    }
    
    

    /**
     * Just reset all learning elements. Note that these are not
     * really gradients as the name of the function implies.
     */
     void LSTDQController::resetGrad() {
	 bVec.clear();
	 lstdqTrace.clear();  
	 AInverse.assign(I); 
	 lstdqSteps = 0;
     }
    
    
    /**
     * The assumption here is that the reward received is for the
     * transition, i.e., the previous state, current state
     */
    void LSTDQController::instantStep(Vector& reward) {

	// Don't know how to do multiple agents here, although it could be done.
	assert(reward.size() == 1);

	// Have to use the previous reward, not the current one.
	lastReward = oldReward;
	oldReward = reward[0];
	instantStepRequired = true;

    }


    /** 
     * Do the work that instantStep() would normally do. This is called from getAction()
     * This assumes rewards are independent of the transition and next state.
     * @param reward the last reward received (between the previous action and current action)
     * @param steps the number of steps done
     */
    void LSTDQController::delayedInstantStep(double reward, int steps) {

	assert(instantStepRequired);

	updateStatistics(reward, steps);

	computeDirection(steps);
    
	// Update the parameters.
	approx->scatter(lstdqWeights, Approximator::PARAMS);
	
	//cout<<"New weights="<<lstdqWeights<<endl;

	
	instantStepRequired = false;

    }


    /**
     * Used in batch mode only. Accumualtes a direction to change without doing a change in params. This is called instead of
     * #instantStep() to implement LSPI. RLAlgs like #GPomdp and #LineSearchAlg use this method.
     */
    void LSTDQController::accumulateGrad(Vector& reward, Observation& newObs) {


	updateStatistics(oldReward, newObs.getSteps());
	oldReward = reward[0];
    }


    /**
     * Update the A matrix (actually, it's inverse) and bVec vector
     * prior to solving the system. 
     */
    void LSTDQController::updateStatistics(double reward, int steps) {

	// Increment the eligibility trace by lambda
	// for real LSTDQ discount=0.0
	lstdqTrace *= discount;
	lstdqTrace += previousPhi;

	//cout<<"Step "<<steps<<endl;
	//cout<<"\tlstdqTrace(lambda="<<discount<<")="<<lstdqTrace<<endl;

	// Add in the observed feature (basis vector) for the "current" state, which 
	// is actually the previous state we observed.
	// previousPhi = phi_{t-1}
	lstdqTD = previousPhi;
	

	// Subtract off the future (i.e., actually the current) observation * discount
	// Note that both phi's depend on the actions that were chosen
	lstdqTD  -= tdDiscount*currentPhi;

	//cout<<"\tlstdqTD="<<lstdqTD<<endl;

	// Okay, we have computed the TD part of the vector. Now do
	// outer product with trace to update A. Although in fact,
	// since what we want in the end is the inverse of A, lets
	// update that directly using the Sherman-Morrison
	// formulation.  If this ever becomes a higher rank update, it
	// uses Sherman-Morrison-Woodbury

	// Alpha implements a moving average of A and b. This stops
	// them blowing up in infinite horizon settings. Be warned
	// that this is probably not the right thing to do while the
	// policy is changing. In that case, decay A and bVec but add
	// in constant new terms.
	double alpha = 1.0 - 1.0/(++lstdqSteps); // First steps value should be 1.0

	if (alpha > 0) AInverse /= alpha; // same as multiplying A by alpha

	lstdqTD *= (1.0 - alpha); // This will weight the new contribution to A by (1.0 - alpha)
	//cout<<"\talpha (lstdqSteps)="<<alpha<<" ("<<lstdqSteps<<")"<<endl;
	//cout<<"\tAInverse="<<AInverse<<endl;
	//cout<<"\tlstdqTD="<<lstdqTD<<endl;

	// If \lambda = 0, then lpiTrace is just the current (i.e.,
	// previous) observation vec embedded in a bunch of 0's
	// depending on the action that was chosen.  Should check that
	// it's not the future observation vec.
	UBlasExtras::updateInverseMatrix(AInverse, lstdqTrace, lstdqTD);
	
	//cout<<"\tA^-1="<<AInverse<<endl;

	bVec *= alpha;
	
	// Update normal (GPOMDP style) gradient
	bVec += reward*(1.0 - alpha)*lstdqTrace;

	//cout<<"\tr="<<reward<<endl;
	//cout<<"\tb="<<bVec<<endl;

    }


    /**
     * Solve the system A^{-1} b = w
     * @param  how many steps did we run to compute grad (not used)
     */
    void LSTDQController::computeDirection(int) {

	// Calculate the new weights lstdqWeights (the w vector)
	axpy_prod(AInverse, bVec, lstdqWeights, true); // true = set lstdqWeights to 0 first

    }


    /**
     * Set the params to the new weights solution.  Unlike the name
     * implies, this is not a step. It's a new absolute solution. The
     * name is reused to fit in with all the other controllers.
     * Now just for batch mode, hence the reset
     */
    void LSTDQController::batchStep() {
   
	// Give the approximator the new solution.
	approx->scatter(lstdqWeights, Approximator::PARAMS);
	resetGrad();

    }


    void LSTDQController::reduce(Vector& v, Approximator::StatsEnum s) {

	switch (s) {
        case Approximator::TRACES: // Deliberate fall through 
        case Approximator::PARAMS: 
	    assert(v.size() == (size_t)getNumParams());
	    v.assign(project(lstdqWeights, range(0, getNumParams())));
	    break;
        case Approximator::GRADS: 
	    approx->reduce(v, s);
	    break;
        default:
	    NYI
		}
    }


    void LSTDQController::scatter(Vector& v, Approximator::StatsEnum s) { 
    
	switch (s) {
        case Approximator::PARAMS: 
	    assert(v.size() == (size_t)getNumParams());
	    project(lstdqWeights, range(0, getNumParams())).assign(v);
	    break;
        case Approximator::TRACES: // Deliberate fall through
        case Approximator::GRADS:
	    approx->scatter(v, s);
	    break;
        default:
	    throw runtime_error("LSTDQController::scatter() unknown type\n");
	}

    }


    /**
     * Takes a local copy of discount and passes it back down to next approx
     */
    void LSTDQController::setDiscount(double discount) {

	this->discount = discount;
	approx->setDiscount(discount);

    }



}
