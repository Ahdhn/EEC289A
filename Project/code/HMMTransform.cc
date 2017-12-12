/**
 * $Id: HMMTransform.cc 127 2007-09-10 17:20:04Z daa $
 */

#include"PGBasics.hh"
#include"HMMTransform.hh"

using namespace std;



namespace libpg {

#ifndef MAXDOUBLE
#define MAXDOUBLE (numeric_limits<double>::max())
#endif

#define RAND_PROB(m) (m) * ( 1. + (random()/(double)RAND_MAX - .5))

HMMTransform::HMMTransform(Controller* controller, 
			   size_t maxHistory, size_t window,
			   size_t numStates, 
			   size_t numObs, 
			   size_t numActions
			   ) : 
    TransformController(controller) {

    // History length should be at least 2, otherwise the program will
    // most probably crash.
    assert(window > 2);
    assert(maxHistory >= window);

    this->maxHistory = maxHistory;
    this->window = window;
    this->numStates = numStates;
    this->numObs = numObs;
    this->numActions = numActions;

    historyLength = 0;
    oHistory.resize(maxHistory); oHistory.clear();
    aHistory.resize(maxHistory); aHistory.clear();
    
    alpha.resize(numStates, maxHistory);    alpha.clear();
     beta.resize(numStates, maxHistory);     beta.clear();
    gamma.resize(numStates, maxHistory);    gamma = scalar_matrix<double> (numStates, maxHistory, 1./numStates);

    xi = new Matrix[numStates];
    for(size_t i=0; i<numStates; i++) {
	xi[i].resize(numStates,maxHistory-1);
	xi[i].clear();
    }

    A = new Matrix[numActions];
    for(size_t a=0; a<numActions; a++) {
	A[a].resize(numStates, numStates);
	for(size_t i=0; i<numStates; i++) {
	    double norm=0.0;
	    for(size_t j=0; j<numStates; j++) {
		A[a](i,j) = RAND_PROB(1./numStates);
		norm += A[a](i,j);
	    }
	    for(size_t j=0; j<numStates; j++)
		A[a](i,j) /= norm;
	}
    }

    B = new Matrix[numActions];
    for(size_t a=0; a<numActions; a++) {
	B[a].resize(numStates, numObs);
	for(size_t i=0; i<numStates; i++) {
	    double norm=0.0;
	    for(size_t j=0; j<numObs; j++) {
		B[a](i,j) = RAND_PROB(1./(numObs));
		norm += B[a](i,j);
	    }
	    for(size_t j=0; j<numObs; j++)
		B[a](i,j) /= norm;
	}
    }

    pi.resize(numStates);
    for(size_t i=0; i<numStates; i++)
	pi(i) = 1./(numStates);

    // HMM vector size must be fixed to input dim of controller.
    // Note that HMM type is Observation.
    hmm.init(controller->getInputDim(), 1, 1, 1);
    hmm.getFeatures().clear();

    //for(int a=0; a < (int) numActions; a++) cerr<<"A["<<a<<"]=\t"<<A[a]<<endl;
    //for(int a=0; a < (int) numActions; a++) cerr<<"B["<<a<<"]=\t"<<B[a]<<endl;

    learn = true;
    prevAction = 0;
}

HMMTransform::~HMMTransform() {
    delete A;
    delete B;
    delete xi;
}

void HMMTransform::getAction(Observation& obs, Vector& action, bool computeGrad) {
    
    // This HMM code can't handle multi-dim obs.
    assert(obs.getFeatures().size1() == 1);

    // Translate histories
    project( oHistory, range(0, maxHistory-1))
        = project( oHistory, range(1, maxHistory));
    project( aHistory, range(0, maxHistory-1))
        = project( aHistory, range(1, maxHistory));
    // Add new observation and action at the end of the queues.
    oHistory(maxHistory-1) = (int)(obs.getFeatures())(0,0);
    aHistory(maxHistory-1) = prevAction;

    // re-estimate the HMM model parameters (when buffer full)
    if (historyLength < maxHistory)
	historyLength++;
    else {
	if ((obs.getSteps()-1)%maxHistory == 0)
	    reEstimate();
	
	estimate();

	/*
	size_t iMax=0; double pMax = gamma(iMax, maxHistory-1);
	for(size_t i = 0; i < numStates; i++)
	    if (gamma(i, maxHistory-1) > pMax) {
		pMax = gamma(i, maxHistory-1);
		iMax = i;
	    }
	gamma.clear();
	gamma(iMax, maxHistory-1) = 1.;
	if (historyLength >= maxHistory) {
	    if (pMax < .5) cerr<<"?"<<endl;
	    else cerr<<iMax<<endl;
	    }
	*/
    }

    // get the final belief state
    // and recompute pi (for next iteration)
    hmm.getFeatures() = project(gamma,
				range(0, numStates),
				range(maxHistory-1, maxHistory)
				);
    //if (historyLength >= maxHistory)
    //cerr<<hmm.features<<endl;

    hmm.setAgent(obs.getAgent());
    hmm.setSteps(obs.getSteps());

    try {
	controller->getAction(hmm, action, computeGrad);
    } catch (std::exception& e) {
	cout<<"hmm: "<<hmm.getFeatures()<<endl;
	throw e;
    }

    // The HMM wants to know what action created what observation, not what
    // observation generates which action.
    prevAction = (size_t)action[0];
    assert(prevAction < numActions);

}

void HMMTransform::reEstimate(void) {
    
    double startP;
    double oldP = -MAXDOUBLE;
    double P=0;
    //int l = 0;
    Vector scale(maxHistory); // initialized by forwardEst()
    
    //cerr<<endl<<endl;
    startP = P = forwardEst(scale, maxHistory); //cerr<<"alpha=\t"<<alpha<<endl;
    
    do {
	//cerr<<endl;
	//cerr<<"HMM loop "<<l++<<": log P="<<P<<endl;
	oldP = P;
	backwardEst(scale, maxHistory); //cerr<<"beta=\t"<<beta<<endl;
	gammaEst(maxHistory); //cerr<<"gamma=\t"<<gamma<<endl;
	xiEst();
	//for(size_t i=0; i<numStates; i++) cerr<<"xi["<<i<<"]=\t"<<xi[i]<<endl;
	maximiseTrans();
	//for(int a=0; a < (int) numActions; a++) cerr<<"A["<<a<<"]=\t"<<A[a]<<endl;
	maximiseDiscreteObs();
	//for(int a=0; a < (int) numActions; a++) cerr<<"B["<<a<<"]=\t"<<B[a]<<endl;
	P = forwardEst(scale,maxHistory); //cerr<<"alpha=\t"<<alpha<<endl;
    }
    while (fabs(P - oldP)/fabs(P) > TERM_PROB_THRESH);
    // Stops if change is less than .1 %
    
    //cerr<<"Restimation took prob from "<<startP<<" --> "<<P<<endl;
}

void HMMTransform::estimate(void) {

    Vector scale(maxHistory); // initialized by forwardEst()

    forwardEst(scale, window); 
    backwardEst(scale, window);
    gammaEst(window);
}

/**
 * Uses log scaling to compute the likelihood that this model fits the
 * observations and rewards
 */
double HMMTransform::forwardEst(Vector &scale, size_t L) {

    double logP;

    size_t i, j; 	/* state indices */
    size_t t;	/* time index */
    size_t T = maxHistory; /* max time */
    
    double sum;	/* partial sum */

    /* 1. Initialization */
    scale(T-L) = 0.0;
    for (i = 0; i < numStates; i++) {
	//alpha(i, 0) = 1.0;             // Assume equally likely to start
	if (T==L)
	    alpha(i, 0) = pi(i) * B[aHistory(0)](i,oHistory(0));
	else
	    alpha(i, T-L) = gamma(i, T-L) * B[aHistory(T-L)](i,oHistory(T-L));
	scale(T-L) += alpha(i, T-L);
    }
    for (i = 0; i <  numStates; i++) {
	alpha(i, T-L) /= scale(T-L); 
    }
    
    /* 2. Induction */
    for (t = T-L; t < T - 1; t++) {
	scale(t + 1) = 0.0;
	for (j=0; j < numStates; j++) {
	    sum = 0.0;
	    for (i = 0; i < numStates; i++) 
		sum += alpha(i, t)*A[aHistory(t+1)](i,j);
	    
	    alpha(j, t + 1) = sum*B[aHistory(t+1)](j,oHistory(t+1));
	    
	    scale(t + 1) += alpha(j, t + 1);
	}
	for (j = 0; j < numStates; j++)
	    alpha(j, t + 1) /= scale(t + 1);
    }
    
    /* 3. Termination */
    logP = 0;
    for (t = T-L; t <  T; t++)
	logP += log(scale(t));
    
    return logP;
}


void HMMTransform::backwardEst(Vector &scale, size_t L) {
    
    size_t i, j;   /* state indices */
    int t;      /* time index */
    double sum;
    size_t T = maxHistory; /* max time */

    /* 1. Initialization */ 
    for (i = 0; i < numStates; i++)
	beta(i, T - 1) = 1.0;
    
    /* 2. Induction */
    for (t = T - 2; t >= (int)(T-L); t--) {
	for (i = 0; i <  numStates; i++) {
	    sum = 0.0;
	    for (j = 0; j <  numStates; j++) 
		sum += A[aHistory(t+1)](i,j)*
		    B[aHistory(t+1)](j,oHistory(t+1))*beta(j,t + 1);
	    
	    beta(i, t) = sum/scale(t);
	}
    }
}

    
/**
 * Gamma is the state occupancy prob at each time step
 */
void HMMTransform::gammaEst(size_t L) {
    
    size_t i;
    size_t t;
    double denominator;
    size_t T = maxHistory;
    
    for (t = T-L; t <  T; t++) {
	denominator = 0.0;
	for (i = 0; i <  numStates ; i++) {
	    gamma(i, t) = alpha(i, t)*beta(i, t);
	    denominator += gamma(i, t);
	}
	
	for (i = 0; i < numStates; i++) 
	    gamma(i, t) = gamma(i, t)/denominator;
    }
    
}


/**
 * xi is the expected transitions from i->j at each time step
 * xi is one shorter than estimationLen since there is one fewer transitions
 * than states occupied.
 */
void HMMTransform::xiEst(void) {

    size_t i, j;
    size_t t;
    double sum;
    size_t T = maxHistory;
    
    for (t = 0; t < T - 1; t++) {
	sum = 0.0;
	for (i = 0; i < numStates   ; i++) {
	    for (j = 0; j < numStates  ; j++) {
		xi[i](j, t) = alpha(i, t)*beta(j, t + 1)*
		    A[aHistory(t+1)](i,j)*B[aHistory(t+1)](j, oHistory(t+1));

		sum += xi[i](j, t);
	    }
	}
	
	// Normalise
	for (i = 0; i < numStates  ; i++) {
	    for (j = 0; j < numStates  ; j++) xi[i](j, t) /= sum;
	}
    }
}




/**
 * Maximise the probability that the observed streams match this models
 * transition parameters.
 */
void HMMTransform::maximiseTrans(void) {

    double denominatorA;
    double numeratorA;
    double newVal;
    size_t i, j;
    size_t T = maxHistory;
    size_t t;

    for(int a = 0; a < (int)numActions; a++)
	// Loop through all transition parameters and update them.
	for (i = 0; i <  numStates; i++) {
	    double sum = 0.;
	    for (j = 0; j <  numStates; j++) {
		denominatorA = 0.0;
		numeratorA = 0.0;
		for (t = 0; t <  T - 1; t++) 
		    if (aHistory(t+1) == a) {
			numeratorA += xi[i](j, t);
			denominatorA += gamma(i, t);
		    }
		
		newVal = MIN_PROB + (1 - MIN_PROB)*numeratorA/denominatorA;
		
		// Look out for situation where denominatorA is 0
		A[a](i,j) = std::isnan(newVal) ? A[a](i,j) : newVal;
		sum += A[a](i,j);
	    }
	    for (j = 0; j <  numStates; j++)
		A[a](i,j) /= sum;
	}
    
}

/**
 * Maximise for discrete symbols
 */
void HMMTransform::maximiseDiscreteObs(void) {

    double denominator;
    Vector numerator(1);
    numerator.resize(this->numObs); // Prevents gcc 4.x warnings
    numerator.clear();

    size_t j;
    size_t T = maxHistory;
    size_t t;
    size_t y;

    for (int a = 0; a < (int)numActions; a++) {
	for (j = 0; j < numStates; j++) {
	    numerator.clear();
	    denominator = 0.0;
	    for (t = 0; t < T; t++)
		if (aHistory(t) == a) {
		    denominator += gamma(j, t);
		    numerator[oHistory(t)] += gamma(j, t);
		}
	    
	    if (denominator != 0) {
		// Update each symbol prob
		double sum = 0.0;
		for (y=0; y < numObs; y++) {
		    B[a](j, y) = MIN_PROB + 
			(1.0 - MIN_PROB)*numerator[y]/denominator;
		    sum += B[a](j,y);
		}
		
		for (y=0; y < numObs; y++)
		    B[a](j,y) /= sum;
	    }
	    else
		for (y=0; y < numObs; y++)
		    B[a](j,y) = 1./numObs;
	}
    }
}
}
