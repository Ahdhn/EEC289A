/**
 * $Id: LookupTable.cc 127 2007-09-10 17:20:04Z daa $
 */

#include"PGBasics.hh"
#include"LookupTable.hh"

using namespace std;
namespace libpg {

/**
 * @param  observations number of discrete observations (not dimensionality) 
 * @param  actions number of discrete actions (not dimensionality) 
 */
LookupTable::LookupTable(int observations, int actions) {

    this->observations = observations;
    this->actions = actions;
    
    params.resize(actions, observations);
    trace.resize(actions, observations);

}


/**
 * Simply return a copy of column associated with obs in output
 * @param  obs observation matrix, in this case always 1x1
 * @param output will be populated with real values action preferences
 */
void LookupTable::doApprox(Observation& obs, Vector& output) { 
    

    int intObs;

    assert(obs.getFeatures().size2() == 1 && obs.getFeatures().size1() == 1);  
    intObs = (int)(obs.getFeatures())(0,obs.getAgent());
    assert(intObs < observations);

    output = column(params, intObs);

}


void LookupTable::feedbackGrad(Observation& obs, Vector& deltas) {
    column(trace, (int)(obs.getFeatures())(0,obs.getAgent())) += deltas;
}


void LookupTable::discountTrace() {
    trace *= discount;
}


void LookupTable::setDiscount(double discount) {
    this->discount = discount;
}


void LookupTable::instantStep(double reward) {

    double multiplier = reward*stepSize;
    params += multiplier*trace;
}


void LookupTable::setStepSize(double stepSize) {
    this->stepSize = stepSize;
}

void LookupTable::resetTrace() {
    trace.clear();
}


void LookupTable::resetParams() {
    params.clear();
}


void LookupTable::randomizeParams(double maxRand) {
    UBlasExtras::randomize(params, maxRand);
}


double LookupTable::getMaxParam() {
    return norm_inf(params);
}


/**
 * Write params to some stream. Usually used for saving paramter values.
 * @param  destination stream
 */
void LookupTable::write(ostream& o) {
    o<<params;
}


/**
 * Write params to some stream. Usually used for saving paramter values.
 * @param  destination stream
 */
void LookupTable::read(istream& o) {
    o>>params;
}


int LookupTable::getNumParams () {
    return params.size1() * params.size2();
}


void LookupTable::scatter(Vector& v, Approximator::StatsEnum s) {
    
    assert(v.size()==(size_t)(params.size1()*params.size2()));
        
    switch (s) {
        case PARAMS:
	    params.clear();
	    UBlasExtras::addScaledVectorToMatrix(1.0, v, params, 0);
	    break;
        case TRACES:
	    params.clear();
	    UBlasExtras::addScaledVectorToMatrix(1.0, v, trace, 0);
	    break;
        default:
	    throw runtime_error("unknown StatsEnum LookupTable::scatter\n");
    }

}


void LookupTable::reduce(Vector& v, Approximator::StatsEnum s) {

    assert(v.size() >= (size_t)params.size1()*params.size2());

    switch (s) {
        case PARAMS:
	    UBlasExtras::addMatrixToVector(params, v, 0);
	    break;
        case TRACES:
	    UBlasExtras::addMatrixToVector(trace, v, 0);
	    break;
        default:
	    throw runtime_error("unknown StatsEnum LookupTable::scatter\n");
    }


}


/**
 * Basic lookup table always has a single dimension input 
 */
int LookupTable::getInputDim() { return 1; }

/**
 * Just the number of actions
 */
int LookupTable::getOutputDim() { return actions; }
}
