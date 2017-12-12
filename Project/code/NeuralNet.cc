/**
 * $Id: NeuralNet.cc 127 2007-09-10 17:20:04Z daa $
 */

#include"PGBasics.hh"
#include"NeuralNet.hh"

using namespace std;
namespace libpg {

/**
 * Create a new linear approximator with given inputs and outputs
 * @param ins input dimensionality
 * @param outs output dimensionality
 */
NeuralNet::NeuralNet(int ins, int outs) {

    dims.resize(2);
    dims[0] = ins;
    dims[1] = outs;

    squash.resize(2);
    squash.clear();

    init();

}


/** Create a new neural network. Creates params, trace variables (to
 * accumulate gradients) activations and error propogation (delta)
 * variables.
 * @param layerDims is a vector that describes the dimension of all
 * layers including input and outputs, in order (el 0 is input layer
 * dim, last el is output dim)
 * @param layerSquash is non-zero if the layer should be
 * squashed. Does not work for input dim.
 */
NeuralNet::NeuralNet(Vector& layerDims, Vector& layerSquash) {

    assert(layerSquash.size() == layerDims.size());
    dims = layerDims;
    squash = layerSquash;
    init();

}


/**
 * After various constructors have figured out layers and sizes, init() is 
 * repsonsible for allocating all the memory
 */
void NeuralNet::init() {

    layers = dims.size() - 1;

    assert(layers > 0); // Need at least one layer
 
    parameters = 0;

    layerParams = new Matrix[layers];
    layerTraces = new Matrix[layers];
    activations = new Vector[layers-1]; // last activation is output
    deltas = new Vector[layers-1];

    // Layer params assume a row of params for each input. Number of
    // columns is number of outputs.
    for (int l=0; l < layers; l++) {
	layerParams[l].resize((int)dims(l), (int)dims(l+1)); // input to output 
	layerTraces[l].resize((int)dims(l), (int)dims(l+1));
     	parameters+=(int)dims[l]*(int)dims[l+1];

	if (l < layers-1) {
	    activations[l].resize((int)dims(l+1)); // same size as layer output
	    deltas[l].resize((int)dims(l+1)); // same size as layer output
	}
    }

}


NeuralNet::~NeuralNet(){
    
    delete[] layerParams;
    delete[] layerTraces;
    delete[] activations;
    delete[] deltas;
    
}


//#define PLSQUASH
#ifdef PLSQUASH
/**
 * piecewise-linear approx to tanh with 5 pieces. about 9 times faster than libc
 * tanh. Not so great for training neural nets, but may be good for
 * implementation
 */
void
NeuralNet::squashVec(Vector& x)
{
    static const double tanh1 = 0.76159415596;
    static const double scale = 0.079245047928; // (tanh(4) - tanh(1)) / (4 - 1)

    Vector::iterator i;
    for (i = x.begin(); i != x.end(); i++) {
        if ((*i) >= -1 && (*i) <= 1)
            (*i) *= tanh1;
        else if ((*i) < -4)
            (*i) = -1;
        else if ((*i) > 4)
            (*i) = 1;
        else if ((*i) < -1)
            (*i) = -tanh1 + ((*i) + 1) * scale;
        else
            (*i) =  tanh1 + ((*i) - 1) * scale;
    }
}

/**
 * Compute v*(dy/dx tanh(x)) = v(1 - y^2)
 */
void NeuralNet::dSquashVec(Vector& y, Vector& v) {

    static const double tanh1 = 0.76159415596;
    static const double scale = 0.079245047928; // (tanh(4) - tanh(1)) / (4 - 1)

    double slope;
    Vector::iterator i;
    for (i = y.begin(); i != y.end(); i++) {
        if ((*i) >= -1 && (*i) <= 1)
            slope = tanh1;
        else if ((*i) < -4 || (*i) > 4)
            slope = 0.;
	else
	    slope = scale;
	v[i.index()] *= slope;
    }

}

#else

/**
 * Compute y = tanh(x)
 */
void NeuralNet::squashVec(Vector& x) {

    for (size_t i = 0; i < x.size(); i++) {
	x[i] = tanh(x[i]);
    }

}

/**
 * Compute v*(dy/dx tanh(x)) = v(1 - y^2)
 */
void NeuralNet::dSquashVec(Vector& x, Vector& v) {

    for (size_t i = 0; i < x.size(); i++) {
	v[i] *= (1.0 - x[i]*x[i]);
    }

}

#endif


/**
 * Feedforward stage of the neural network. If you don't know how
 * Neural nets work, don't even bother trying to understand this bit
 * of the code.
 * @param  input observation
 * @param  output approximation destination matrix
 */
void NeuralNet::doApprox(Observation& obs, Vector& output) {

    assert(obs.getFeatures().size1() == (unsigned int)dims(0));
    // For efficiency (both space and time) we treat the single layer
    // case and multi-layer cases differently, using the obs as the
    // first activiation and output as the final output activation.

    //cout<<"Obs="<<obs.getFeatures()<<endl;
    //cout<<layerParams[0];

    //Linear network
    if (layers < 2) {
	axpy_prod(column(obs.getFeatures(), obs.getAgent()), layerParams[0], output, true);
	if ((int)squash[0])  { squashVec(output); }
		
    } 
    else {
	int l=0;

	// First layer
	axpy_prod(column(obs.getFeatures(), obs.getAgent()), layerParams[l],  activations[l], true);
	// Squash output of first layer? A little confusing since we
	// must check squash[l+1] to see if we squash activations[l]
	if ((int)squash[l+1]) squashVec(activations[l]);
	l++;
	
	// Middle layers
	for (; l < layers-1; l++) {
	    axpy_prod(activations[l-1], layerParams[l], activations[l], true);
	    if ((int)squash[l+1]) squashVec(activations[l]);
	}
	
	// Output layer
	axpy_prod(activations[l-1], layerParams[l], output, true);
	if ((int)squash[l+1]) squashVec(output);

    }

}


/**
 * Do error back propogation, but on the outputDeltas provided rather
 * than on any error
 * @param  input observation (again)
 * @param  output deltas, i.e., the gradient of the softmax function.
 * SIDEEFFECT: adds gradients to the trace variables in this class
 */
void NeuralNet::feedbackGrad(Observation& obs, Vector& outputDeltas) {

    // Single layer case
    if (layers < 2) {
	// rank one update = outer produc of two vectors
	// Compute grad for only layer
	noalias(layerTraces[0]) += outer_prod(column(obs.getFeatures(), obs.getAgent()), outputDeltas);
    }
    else {
	int l = layers - 1;

	// Compute grad for last layer
	noalias(layerTraces[l]) += outer_prod(activations[l-1], outputDeltas);
	// Back prop deltas through last layer
	axpy_prod(layerParams[l], outputDeltas, deltas[l-1], true);
	
	// Feedback through squashing
	// activations[0] is OUTPUT of first layer
	if ((int)squash[l]) dSquashVec(activations[l-1], deltas[l-1]);
	l--;

	// repeat for all hidden layers
	for (; l > 0; l--) {
	    // Compute grad for current layer
	    noalias(layerTraces[l]) += outer_prod(activations[l-1], deltas[l]);
	    // Back prop deltas through current layer
	    axpy_prod(layerParams[l], deltas[l], deltas[l-1], true);
	    // Feedback through squashing
	    if ((int)squash[l]) dSquashVec(activations[l-1], deltas[l-1]);
	}

	// Compute grad for input layer
	noalias(layerTraces[l]) += outer_prod(column(obs.getFeatures(), obs.getAgent()), deltas[l]);
	
    }

}


/**
 * Normal discounting, but have to go through all the para matricies
 */
void NeuralNet::discountTrace() {
    for (int l=0; l < layers; l++) {
	layerTraces[l] *= discount;
    }

}


/**
 * Set the discount factor for discounting the traces
 */
void NeuralNet::setDiscount(double discount){ 
    this->discount = discount;
}


/**
 * Step in the direction of the trace*reward. Note that we
 * pre-axpy_prodiply the reward and step size which saves a lot of work
 * over a naive implementation.
 * @param  reward for this step.
 */
void NeuralNet::instantStep(double reward){
 
    double multiplier = reward*stepSize;
    for (int l=0; l < layers; l++) {
	noalias(layerParams[l]) += multiplier*layerTraces[l];
    }
}



/**
 * Set the step size for OLPomdp. This can evolve over time if
 * controlled by OLPomdp @param  new step size
 */
void NeuralNet::setStepSize(double stepSize){ 
    this->stepSize = stepSize;
}


/**
 * Set all traces to 0, e.g., starting learning from scratch, or the
 * end of an episode.
 */
void NeuralNet::resetTrace() { 
    for (int l=0; l < layers; l++) {
	layerTraces[l].clear();
    }

}


/**
 * Set all params to 0. equivalent to a uniform policy.  This is a
 * really bad idea for multi-layer perceptrons. For MLPs, use small
 * random values instead, and pray that I've written that function by
 * the time you want to use it.
 */
void NeuralNet::resetParams() { 
    for (int l=0; l < layers; l++) {
	layerParams[l].clear();
    }
}


/**
 * Randomize the paramters, usually to a small value. Important for
 * multi-layer neural nets
 * @param  abs max value of random number.
 */
void NeuralNet::randomizeParams(double maxRand) { 

  for (int l=0; l < layers; l++) {
	UBlasExtras::randomize(layerParams[l], maxRand);
    }
 
}


/**
 * @return  maximum absolute parameter value
 */
double NeuralNet::getMaxParam() {
    
    double maxi = 0.0;
    for (int l=0; l < layers; l++) {
	maxi = std::max(maxi, (double)norm_inf(layerParams[l]));
    }

    return maxi;
}


/**
 * Write params to some stream. Usually used for saving paramter values.
 * @param  destination stream
 */
void NeuralNet::write(ostream& o) {

    for (int l=0; l < layers; l++) {
	o<<layerParams[l]<<endl;
    }

}


/**
 * Write params to some stream. Usually used for saving paramter values.
 * @param  destination stream
 */
void NeuralNet::read(istream& o) {

    for (int l=0; l < layers; l++) {
	o>>layerParams[l];
    }

}


int NeuralNet::getNumParams() {
    return parameters;
}


int NeuralNet::getInputDim() {
    return (int)dims[0];
}


int NeuralNet::getOutputDim() {
    return (int)dims[layers];
}


/**
 * @param s whether to collect params or  traces. GRADS is invalid
 */
void NeuralNet::reduce(Vector& v, StatsEnum s) {

    assert(v.size() >= (size_t)parameters);

    int p = 0;

    for (int l=0; l < layers; l++) {
	switch (s) {
	    case PARAMS: 
		UBlasExtras::addMatrixToVector(layerParams[l], v, p);
		break;
	    case TRACES:
		UBlasExtras::addMatrixToVector(layerTraces[l], v, p);
		break;
	    default:
		throw runtime_error("NeuralNet::reduce() unknown type\n");
	}
	p += layerParams[l].size1()*layerTraces[l].size2();
    }    
    assert(p == parameters);
}


/**
 * @param s what to set, can't do GRADS here.
 */
void NeuralNet::scatter(Vector& v, StatsEnum s) {

    assert(v.size()==(size_t)parameters);

    int p = 0;

    for (int l=0; l < layers; l++) {
	
	switch (s) {
	    case PARAMS:
		layerParams[l].clear();
		UBlasExtras::addScaledVectorToMatrix(1.0, v, layerParams[l], p);
		break;
	    case TRACES:
		layerTraces[l].clear();
		UBlasExtras::addScaledVectorToMatrix(1.0, v, layerTraces[l], p);
	        break;
	    default:
		throw runtime_error("NeuralNet::reduce() unknown type\n");
	}
	    
	p += layerParams[l].size1()*layerParams[l].size2();
    }    
    assert(p == parameters);

}
}
