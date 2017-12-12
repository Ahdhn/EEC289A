/**
 * $Id: LineSearchAlg.cc 127 2007-09-10 17:20:04Z daa $
 */


#include"PGBasics.hh"
#include"GPomdp.hh"
#include"LineSearchAlg.hh"

using namespace std;

namespace libpg {

LineSearchAlg::LineSearchAlg(Controller* controller, 
			     Simulator* simulator, 
			     double discount, 
			     double stepSize, 
			     int lineSearchSteps,  
			     double tolerance, 
			     double backOff,
			     double exponentBase,
			     double downhilllThresh,
			     int maxTries) : GPomdp(controller, simulator, discount, stepSize) {

    lineSearchInit(lineSearchSteps, tolerance, backOff, exponentBase, downhilllThresh, maxTries);

}


void LineSearchAlg::lineSearchInit( int lineSearchSteps,  
				    double tolerance, 
				    double backOff,
				    double exponentBase,
				    double downhilllThresh,
				    int maxTries) {
    
    // Additional parameters needed to control the line search.
    this->lineSearchSteps = lineSearchSteps;
    this->tolerance = tolerance;
    this->backOff = backOff;
    this->exponentBase = exponentBase;
    this->downhilllThresh = downhilllThresh;
    this->maxTries = maxTries;
    
    stepSizeBackup = stepSize;
    
    origParams.resize(controller->getNumParams());
    params.resize(controller->getNumParams());
    grads.resize(controller->getNumParams());
    testGrads.resize(controller->getNumParams());
    
    seed=0;
}



/**
 * @return  totalSteps counts the number of steps we performed evaluating different step sizes.
 */
double LineSearchAlg::lineSearch(int& totalSteps, double lastEpochVal, Vector& totalRewards) {

    double currentStep = stepSize;
    double bestStep = 0.0;
    double currentReward = lastEpochVal;
    double bestReward = lastEpochVal;
    double lastIP = 0; // Last inner product
    double factor = exponentBase;
    int its = 1;
    bool backward = false;

    // To begin, lets get a copy of the current parameters. This is expensive. Don't do it too often.
    origParams.clear();
    controller->reduce(origParams, Approximator::PARAMS);
    grads.clear();
    controller->reduce(grads, Approximator::GRADS);


    // Immediate exit if grad is too small.
    if (norm_2(grads) < tolerance) {
	maxSteps = 1;
	cout<<"\tLS: R="<<bestReward<<" : Terminating because gradient<tolerance\n";
	return bestReward;
    }
    
    lastIP = incrementLineSearch(currentStep, totalRewards, totalSteps, currentReward);
        
    if (lastIP < -tolerance) {
	// We have gone backwards. Reset step size and program to
	// exponentially reduce step size.
	cout<<"\tLS: its=0 R="<<currentReward<<" stepSize="<<currentStep<<" IP="<<lastIP<<" : heading backward\n";
	backward = true;
	currentStep = stepSize/factor;
	lastIP = incrementLineSearch(currentStep, totalRewards, totalSteps, currentReward);
    }
    
    
    if (fabs(lastIP) < tolerance) {
	// We've either magically found the local optimum, or more
	// likely we've done something stupid and gone to an
	// extremum. Try resetting the step size to it's inital
	// value. Assumes this is a small value.
	cout<<"\tLS: its=0 R="
	    <<currentReward
	    <<" stepSize="
	    <<currentStep
	    <<" IP="
	    <<lastIP
	    <<" : reset after immediate extremum\n";
	stepSize = stepSizeBackup;
	currentStep = stepSize;
	lastIP = incrementLineSearch(currentStep, totalRewards, totalSteps, currentReward);
    }

    cout<<"\tLS: its=1 R="<<currentReward<<" stepSize="<<currentStep<<" IP="<<lastIP;

    while (its < maxTries && 
	   ((!backward && lastIP > tolerance) || (backward && lastIP < -tolerance))) {
	if (currentReward > bestReward) {
	    bestReward = currentReward;
	    bestStep = currentStep;
	    cout<<" *";
	}
	cout<<endl;
	currentStep *= (backward?1.0/factor:factor);
	lastIP = incrementLineSearch(currentStep, totalRewards, totalSteps, currentReward);
	cout<<"\tLS: its="<<++its<<" R="<<currentReward<<" stepSize="<<currentStep<<" IP="<<lastIP;       
    }
    cout<<endl;

    if  ((backward && lastIP > tolerance) || (!backward && lastIP < -tolerance)) {
      // We stopped because ip changed sign, within tolerance. 
      // Backoff step to before ip changed sign. 
      // Send out nre parameters.
      currentStep /= (backward?1.0/factor:factor);
      if (currentStep == bestStep) currentReward =  bestReward;
      params.assign(origParams + currentStep*grads);
      controller->scatter(params, Approximator::PARAMS);  
    }

    // Finished line search. Do sanity check on final reward value.
    if (bestReward > currentReward*(downhilllThresh+1)) {
	// Not allowed to make no progress. But be conservative if all else fails.
	if (bestStep == 0) bestStep = stepSizeBackup; 
	params.assign(origParams + bestStep*grads);
	controller->scatter(params, Approximator::PARAMS);   
	currentReward = bestReward;
	currentStep = bestStep;
    }
    
    stepSize = currentStep/backOff;

    cout<<"\tLS: its="<<its++<<" R="<<currentReward<<" stepSize="<<currentStep<<" IP="<<lastIP<<" : Final"<<endl;
    return bestReward;


}



double LineSearchAlg::incrementLineSearch(double currentStep, Vector& totalRewards, int& totalSteps, double& currentReward) {

    // Try the current stepSize   
    params.assign(origParams + currentStep*grads);
    controller->scatter(params, Approximator::PARAMS);   
    controller->resetTrace(); // Normally the batch estimates continue with previous trace.
    srandom(seed);
    currentReward = doSteps(totalRewards, lineSearchSteps, true);
    controller->computeDirection(lineSearchSteps);
    totalSteps += lineSearchSteps;
    testGrads.clear();
    controller->reduce(testGrads, Approximator::GRADS);
    // If inner prod between test grad and start grad is negative, we
    // are going downhill
    return inner_prod(testGrads, grads);
    
}
 

double LineSearchAlg::learnCore(int stepsPerEpoch, int& totalSteps) {

    double lastEpochVal;
    Vector totalRewards(simulator->getRewardDim());

    // Reset random number generator to help line search
    srandom(seed);
    
    // Estimate grad
    lastEpochVal = doSteps(totalRewards, stepsPerEpoch, true);
    totalSteps += stepsPerEpoch;
    
    // Estimate search direction
    controller->computeDirection(stepsPerEpoch);
    
    printPerformanceInfo(false, totalSteps,  lastEpochVal,  controller->getMaxParam(), 0);
    lastEpochVal = lineSearch(totalSteps, lastEpochVal, totalRewards);

    if (useAutoBaseline) baseline.assign(totalRewards);

    return lastEpochVal;

}

}
