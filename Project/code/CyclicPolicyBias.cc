/**
 * $Id: CyclicPolicyBias.cc 127 2007-09-10 17:20:04Z daa $
 */

#include"PGBasics.hh"
#include"CyclicPolicyBias.hh"

using namespace std;
namespace libpg {

/**
 * See header for purpose of this approximator. Always sits on top of
 * another controller that has parameters.
 */
CyclicPolicyBias::CyclicPolicyBias(Approximator* approx, int cycleTime, double bias) {

    this->approx = approx;
    this->bias = bias;
    this->cycleTime = cycleTime;

    if (cycleTime%approx->getOutputDim() != 0) {
	cout<<"!! CyclicPolicyBias requires cycleTime="
	    <<cycleTime
	    <<" have outputs="
	    <<approx->getOutputDim()
	    <<" as a factor\n";
	exit(EXIT_FAILURE);
    }
    
    stepsPerControl = cycleTime/approx->getOutputDim();

}


/**
 * Adds constant to normal approximator
 */
void CyclicPolicyBias::doApprox(Observation& obs, Vector& output) {

    approx->doApprox(obs, output);
    
    output[(obs.getSteps()%cycleTime) / stepsPerControl] += bias;
}


/**
 * Everything else is just a pass through to the normal approximator
 */
void CyclicPolicyBias::feedbackGrad(Observation& obs, Vector& deltas) {
    approx->feedbackGrad(obs, deltas);
}

void CyclicPolicyBias::discountTrace() {
    approx->discountTrace();
}

void CyclicPolicyBias::setDiscount(double discount) {
    approx->setDiscount(discount);
}

void CyclicPolicyBias::instantStep(double reward) {
    approx->instantStep(reward);
}

void CyclicPolicyBias::setStepSize(double stepSize) {
    approx->setStepSize(stepSize);
}

void CyclicPolicyBias::resetTrace() {
    approx->resetTrace();
}

void CyclicPolicyBias::resetParams() {
    approx->resetParams();
}

void CyclicPolicyBias::randomizeParams(double maxRand) {
    approx->randomizeParams(maxRand); 
}

double CyclicPolicyBias::getMaxParam() {
    return approx->getMaxParam();
}

void CyclicPolicyBias::write(ostream& o) {
    approx->write(o);
}

void CyclicPolicyBias::read(istream& o) {
    approx->read(o);
}

void CyclicPolicyBias::batchStep() {
    approx->batchStep();
}

void CyclicPolicyBias::accumulateGrad(double reward, Observation& newObs) {
    approx->accumulateGrad(reward, newObs);
}

void CyclicPolicyBias::resetGrad() {
    approx->resetGrad();
}

void CyclicPolicyBias::computeDirection(int steps) {
    approx->computeDirection(steps);
}


int CyclicPolicyBias::getNumParams() {
    return approx->getNumParams();
}


int CyclicPolicyBias::getInputDim() {
    return approx->getInputDim();
}


int CyclicPolicyBias::getOutputDim() {
    return approx->getOutputDim();
}


void CyclicPolicyBias::reduce(Vector&v, Approximator::StatsEnum s) {
    approx->reduce(v, s);
}


void CyclicPolicyBias::scatter(Vector&v, Approximator::StatsEnum s) {
    approx->scatter(v, s);
}
}
