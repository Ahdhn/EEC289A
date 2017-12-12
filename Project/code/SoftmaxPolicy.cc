#include"PGBasics.hh"
#include <float.h>
#include"SoftmaxPolicy.hh"
#include"Sampler.hh"

using namespace std;

namespace libpg {
	
    SoftmaxPolicy::SoftmaxPolicy(double constTemp) {
	assert(constTemp > 0.0);
	temp = constTemp;
	tempFunc = NULL;
    }

    SoftmaxPolicy::SoftmaxPolicy(Temperature* tempFunc) {
	this->tempFunc = tempFunc;
    }

    SoftmaxPolicy::~SoftmaxPolicy() {
	if (tempFunc != NULL)
	    delete tempFunc;
    }

    void SoftmaxPolicy::getAction(Vector dist, Vector& action, Observation& obs) {
	int sampledAct;
	double oneNorm = 0.0;

	// Computes new temperature according to decaying function
	if (tempFunc != NULL)
	    temp = max(tempFunc->getValue(obs.getSteps()), DBL_MIN);
 
	// Computes Gibbs distribution
	for (Vector::iterator i = dist.begin(); i != dist.end(); i++) {
	    *i = min(exp(min(*i/temp, DBL_MAX)), DBL_MAX);
	    oneNorm = min(oneNorm + *i, DBL_MAX);
	}

	// Pick an action
	sampledAct = Sampler::discrete(dist, oneNorm);
	assert(sampledAct < (int)dist.size());
	action[0] = sampledAct;
    }
}
