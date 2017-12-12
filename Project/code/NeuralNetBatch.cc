/**
 * $Id: NeuralNetBatch.cc 127 2007-09-10 17:20:04Z daa $
 */

#include"PGBasics.hh"
#include"NeuralNet.hh"
#include"NeuralNetBatch.hh"

namespace libpg {
/**
 * A version of neural network that supports aggregated gradients
 */


/**
 * Create a new linear approximator with given inputs and outputs
 * @param ins input dimensionality
 * @param outs output dimensionality
 */
NeuralNetBatch::NeuralNetBatch(int ins, int outs) : NeuralNet(ins, outs) {

    grads = new Matrix[layers+1];

    for (int l=0; l < layers; l++) {
	grads[l].resize((int)dims[l], (int)dims[l+1]); 
    }

}


/**
 * Same as neural network call, but additionaly resizes accumulated
 * gradient matrices
 */
NeuralNetBatch::NeuralNetBatch(Vector& dims, Vector& squash) : NeuralNet(dims, squash) {

    grads = new Matrix[layers+1];

    for (int l=0; l < layers; l++) {
	grads[l].resize((int)dims[l], (int)dims[l+1]); 
    }
}
 

   
void NeuralNetBatch::batchStep() {
    
    for (int l=0; l < layers; l++) {
	noalias(layerParams[l]) += stepSize*grads[l];
    }
}



void NeuralNetBatch::computeDirection(int steps) {
    
    for (int l=0; l < layers; l++) {
	grads[l]/=steps;
    }
}

    
void NeuralNetBatch::accumulateGrad(double reward, Observation& newObs) {
    for (int l=0; l < layers; l++) {
	noalias(grads[l]) += reward*layerTraces[l];
    }
}


void NeuralNetBatch::resetGrad() {
    for (int l=0; l < layers; l++) {
	grads[l].clear();
    }
}
    
/**
 * @param s whether to collect params or grads or traces. Only grads
 * done in this routine. Others passed to NeuralNet
 */
void NeuralNetBatch::reduce(Vector& v, StatsEnum s) {

    assert(v.size() >= (size_t)parameters);

    int p = 0;

    if (s == Approximator::GRADS) {
	// Doing grads
	for (int l=0; l < layers; l++) {
	    UBlasExtras::addMatrixToVector(grads[l], v, p);
	    p += layerParams[l].size1()*layerParams[l].size2();
	}    
	assert(p == parameters);
    }
    // Neural net knows how to do other types
    else NeuralNet::reduce(v, s);

}


void NeuralNetBatch::scatter(Vector& v, StatsEnum s) {

    assert(v.size()==(size_t)parameters);

    int p = 0;

    if (s==Approximator::GRADS) {
	for (int l=0; l < layers; l++) {
	    grads[l].clear();
	    UBlasExtras::addScaledVectorToMatrix(1.0, v, grads[l], p);
	    p += layerParams[l].size1()*layerParams[l].size2();
	}    
	assert(p == parameters);
    }
    else NeuralNet::scatter(v, s);

}
}
