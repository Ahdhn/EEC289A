/*
 *  MultidimensionalLookupTable.cc
 *  
 *
 *  Created by Owen Thomas on 20/06/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */
#include "PGBasics.hh"
#include "LookupTable.hh"
#include "LookupTableBatch.hh"
#include "MultidimensionalLookupTable.hh"

namespace libpg {

/**
 *Integer power function returns y raised to the x.
 */
int power (int y, unsigned int x) {
    int result = 1;
    while(x-- > 0) result *= y;
    return result;
}


/**
 * Construct a FiniteMemoryLookup table.
 *
 * @param observations, the range of observation values.
 * @param outputs, the number of outputs (actions) to map the observations to.
 * @param memoryLength, the number of past observations being remembered.
 */
MultidimensionalLookupTable::MultidimensionalLookupTable (int observations, 
							  int outputs, 
							  int maxHistory) 
{
    
    assert(maxHistory > 0);
    
    this->observations = observations;
    actions = outputs;
    this->maxHistory = maxHistory;
    
    params.resize (outputs,  power (observations, maxHistory)); params.clear();
    trace.resize (outputs, power (observations, maxHistory)); trace.clear();
    grads.resize (outputs, power (observations, maxHistory)); grads.clear();
}


void MultidimensionalLookupTable::doApprox (Observation& obs, Vector& output) {
    
    assert(obs.getFeatures().size1() == (size_t)maxHistory);
    assert(obs.getFeatures().size2() == 1);
    
    int index = 0;
    for (int i = maxHistory-1; i >= 0; i--) {
	index = index*observations + (int)obs.getFeatures()(i, obs.getAgent());
    }
    assert (index < power(observations, maxHistory));
    
    output = column (params, index);
    
}

void MultidimensionalLookupTable::feedbackGrad (Observation& obs, Vector& deltas) {
    
    int index = 0;
    for (int i = maxHistory-1; i >= 0; i--) {
	index = index*observations + (int)obs.getFeatures()(i, obs.getAgent());
    }
    assert (index < power(observations, maxHistory));
    column(trace, index) += deltas;
}



/**
 * An input or each element of history
 */
int MultidimensionalLookupTable::getInputDim() { return maxHistory; }
}
