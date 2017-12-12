/**
 * $Id: LookupTableBatch.cc 127 2007-09-10 17:20:04Z daa $
 */

#include"PGBasics.hh"
#include"LookupTable.hh"
#include"LookupTableBatch.hh"

namespace libpg {

/**
 * A version of the lookup table that supports aggregated gradients
 */


/**
 * Same as neural network call, but additionaly resizes accumulated
 * gradient matrices
 */
LookupTableBatch::LookupTableBatch(int observations, int outputs) : LookupTable(observations, outputs) {

    grads.resize(outputs, observations); 
}
 
   
void LookupTableBatch::batchStep() {
    noalias(params) += stepSize*grads;
}

void LookupTableBatch::computeDirection(int steps) {
    grads/=steps;
}
    
void LookupTableBatch::accumulateGrad(double reward, Observation& newObs) {
    noalias(grads) += reward*trace;
}


void LookupTableBatch::resetGrad() {
    grads.clear();
}


void LookupTableBatch::scatter(Vector& v, Approximator::StatsEnum s) {
    
    assert(v.size()==(size_t)(params.size1()*params.size2()));
        
    switch (s) {
        case GRADS:
	    params.clear();
	    UBlasExtras::addScaledVectorToMatrix(1.0, v, grads, 0);
	    break;
        default:
	    LookupTableBatch::scatter(v, s);
    }

}


void LookupTableBatch::reduce(Vector& v, Approximator::StatsEnum s) {

    assert(v.size() >= (size_t)params.size1()*params.size2());

    switch (s) {
        case GRADS:
	    UBlasExtras::addMatrixToVector(grads, v, 0);
	    break;
        default:
	    LookupTable::reduce(v, s);
    }

}
}
