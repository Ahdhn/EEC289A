
#include <float.h>
#include <math.h>

#include"PGBasics.hh"
#include"eGreedyPolicy.hh"


using namespace std;

namespace libpg {

    eGreedyPolicy::eGreedyPolicy(double e) {
	assert((e >= 0.0) && (e <= 1.0));
	this->e = e;
	epsilonFunc = NULL;
    }


    eGreedyPolicy::eGreedyPolicy(EpsilonFunction* epsilonFunc) {
	this->epsilonFunc = epsilonFunc;
    }


    eGreedyPolicy::~eGreedyPolicy() {
	if (epsilonFunc != NULL)
	    delete epsilonFunc;
    }


    void eGreedyPolicy::getAction(Vector dist, Vector& action, Observation& obs) {
	if (epsilonFunc != NULL)
	    e = epsilonFunc->getValue(obs.getSteps());
	
	if (random()/(double)RAND_MAX >= e) {
	    // Performs greedly
	    double maxValue = dist[0];
	    action[0] = 0;
	    for (Vector::iterator i = dist.begin(); i != dist.end(); i++)
		if (*i > maxValue) {
		    maxValue = *i;
		    action[0] = i.index();
		}
	}
	else
	    // Performs randomly
	    action[0] = floor((random()/(double)RAND_MAX) * dist.size());

	assert(action[0] < (int)dist.size());
    }
}
