/**
 * $Id: NACTransform.cc 127 2007-09-10 17:20:04Z daa $
 */

#include"PGBasics.hh"
#include"NeuralNet.hh"
#include"NeuralNetBatch.hh"
#include"NACTransform.hh"

using namespace std;
namespace libpg {

/**
 * Constructor. Set up additional vectors and matrices needed to
 * implement natural actor critic.
 */
NACTransform::NACTransform(Controller* controller, double tdDiscount) 
    : TransformController(controller) {

 
    this->tdDiscount = tdDiscount;

    // We turn int ovbservations into unit vectors for first term
    // Second term comes from number of parameters.
    nacDimension = controller->getInputDim() + controller->getNumParams();

    I.resize(nacDimension, nacDimension);
    I.clear();
    //A.resize(nacDimension, nacDimension);
    AInverse.resize(nacDimension, nacDimension);
    AInverse.clear();
 
    nacTD.resize(nacDimension);
    nacTrace.resize(nacDimension);
    nacTrace.clear();
    nacGrad.resize(nacDimension);
    nacGrad.clear();
    nacDirection.resize(nacDimension);
    nacDirection.clear();
    nacParamUpdate.resize(controller->getNumParams());
    nacParamUpdate.clear();
    I.assign(identity_matrix<double>(nacDimension));
         
    instantStepRequired = false;


}

 
/**
 * NAC delays it's discounting till the end of a updateStatistics()
 */
void NACTransform::discountTrace() {

}


/**
 * Just a wrapper that intercepts the new observation and completes
 * the instant step from the last round if necessary.
 */
void NACTransform::getAction(Observation& obs, Vector& action, bool computeGrad) {

    assert(obs.getFeatures().size1() == (size_t)getInputDim());
    
    if (instantStepRequired) delayedInstantStep(lastReward, obs);

    controller->resetTrace();    
    controller->getAction(obs, action, computeGrad);

    if (computeGrad) {

	// Grab the most recent trace and add to NAC's trace. Will discount later.
	controller->reduce(nacTrace, Approximator::TRACES);
	
	// Add in the observed feature (basis vector) for the current
	// state into the bottom part of the nacTrace vector.
	// comment this line out if want to use natrual gradient
	project(nacTrace, range(controller->getNumParams(), nacDimension)) += column(obs.getFeatures(), obs.getAgent());
	
	// Make a copy of the last observation.
	lastObs = obs;
    }
}


void NACTransform::resetTrace() {
    nacTrace.clear();
    AInverse.assign(I); 
}


void NACTransform::resetGrad() {
    nacGrad.clear();
    nacTrace.clear();  
    AInverse.assign(I); 
 
}



/**
 * In the case of NAC this is a bit of a hack. The update requires the next 
 * observation to work, so we'll store the fact that an update is required
 * and do it in the next feedback() phase. This effectively means you end
 * up with steps-1 updates, but oh well...
 */
void NACTransform::instantStep(Vector& reward) {

    if (reward.size() > 1) throw runtime_error("NACTransform::instantStep() requires scalar reward\n");
    instantStepRequired = true;
    lastReward = reward[0];

}

void NACTransform::delayedInstantStep(double reward, Observation& newObs) {


    assert(instantStepRequired);
    assert(newObs.getFeatures().size1() == (size_t)getInputDim());

    nacGrad.clear();
    
    updateStatistics(reward, newObs);

    // Divide gradient by A, effectllively applying second order statistics
    axpy_prod(AInverse, nacGrad, nacDirection, true); 
    
    // Update the parameters.
    batchStep();

    instantStepRequired = false;

}


void NACTransform::accumulateGrad(double reward, Observation& newObs) {

    updateStatistics(reward, newObs);
}


void NACTransform::updateStatistics(double reward, Observation& newObs) {
	
    assert(newObs.getFeatures().size1() == (size_t)getInputDim());

    // Construct the nacTD vector 
    nacTD.clear();

    controller->reduce(nacTD, Approximator::TRACES);


    // Add in the observed feature (basis vector) for the current
    // state into the bottom part of the nacTD vector.
    // comment this line out if want to use natural gradient
    project(nacTD, range(controller->getNumParams(), nacDimension)) = 
	column(lastObs.getFeatures(), lastObs.getAgent());

	    
    // Subtract off the future observation * discount
    // comment this line out if want to use natural gradient
    project(nacTD, range(controller->getNumParams(), nacDimension)) -= 
	tdDiscount*column(newObs.getFeatures(), newObs.getAgent());


    // Okay, we have computed the TD part of the vector. Now do outer
    // product with trace to update A Although in fact, since what we
    // want in the end is the inverse of A, lets update that directly.

    double alpha = 1.0 - 1.0/(lastObs.getSteps()+1);
    assert(alpha > 0.0);
    
    //A += (1.0-alpha)*(outer_prod(nacTD, nacTD)+1e-2*I);//natural gradient ggT
    //A += (1.0-alpha)*(outer_prod(nacTrace, nacTD)+1e-2*I);
    //A += outer_prod(nacTrace, nacTD); //it is easily to blow up at some stage
    //A += outer_prod(nacTrace, nacTD)+1e-2*identity_matrix<double>(nacDimension);
    //cout<<"A"<<A<<endl;


    // Update inverse of A
    AInverse /= alpha; // same as multiplying A by alpha)
    nacTD *= (1.0 - alpha);
    UBlasExtras::updateInverseMatrix(AInverse, nacTrace, nacTD);
	
	
    // Update normal (GPOMDP style) gradient
    nacGrad += reward*nacTrace;

    // Update the trace out of order (normal discountTrace command ignored).
    nacTrace *= discount;
}


/**
 * Normalise the gradient and multiply by the inverse of the
 * A matrix, doing the magic step of natural actor-critic in one
 * go.
 * Could use Sherman-Morrison:
 * (A + u^T v)^-1 = A^-1 - (A^-1 u)^T (v A^-1) / (1 + v A^-1 u)
 * @param  how many steps did we run to compute grad
 */
void NACTransform::computeDirection(int steps) {
    
    nacGrad/=steps;

    axpy_prod(AInverse, nacGrad, nacDirection, true); 
    //cout<<"direction="<<nacDirection<<endl;

}


void NACTransform::batchStep() {
    
    
    nacParamUpdate.assign(project(nacDirection, range(0, controller->getNumParams())));
    // Loop over each layer of the neural network
    controller->scatter(nacParamUpdate, Approximator::GRADS);
    controller->batchStep();
    
}


void NACTransform::reduce(Vector& v, Approximator::StatsEnum s) {

    switch (s) {
        case Approximator::TRACES: // Deliberate fall through 
        case Approximator::PARAMS: 
	    controller->reduce(v, s);
	    break;
        case Approximator::GRADS: 
	    assert(v.size() == (size_t)getNumParams());
	    v.assign(project(nacDirection, range(0, getNumParams())));
        default:
	    NYI
    }
}


void NACTransform::scatter(Vector& v, Approximator::StatsEnum s) { 
    
    switch (s) {
        case Approximator::PARAMS: // deliberate fall through
        case Approximator::TRACES:
	    controller->scatter(v, s);
	    break;
        case Approximator::GRADS:
	    assert(v.size() == (size_t)getNumParams());
	    project(nacDirection, range(0, getNumParams())).assign(v);
        default:
	    throw runtime_error("NACTransform::scatter() unknown type\n");
    }

}


/**
 * Takes a local copy of discount and passes it back down to next controller
 */
void NACTransform::setDiscount(double discount) {

    this->discount = discount;
    controller->setDiscount(discount);

}
}
