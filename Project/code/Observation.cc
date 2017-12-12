#include "PGBasics.hh"

namespace libpg {

    Observation::Observation(int rows, int cols, int agents, int actions) {

	init(rows, cols, agents, actions);

    }


    /**
     * Copy constructor
     */
    Observation::Observation(Observation& obs) {

	*this = obs;
    }


    Observation& Observation::operator=(Observation& rhs) {
	if (&rhs != this) {
	    features = rhs.getFeatures();
	    eligible = rhs.getEligible();
	    agent = rhs.getAgent();
	    steps = rhs.getSteps();
	}
	return *this;
    }


    void Observation::init(int rows, int cols, int agents, int actions) {

	// Arbitrary matrix
	features.resize(rows, cols);
	
	// In the single agent case the eligible vector says which actions
	// are eligible. In the multi-agent case it says which agents are
	// eligible.
	eligible.resize((agents==1)?actions:agents);
	eligible.clear();
  
	// This is initialised to all ones so that the default behaviour is
	// to run all agents. This means that if you don't know anything
	// about the eligible vector you should still get the expected
	// behaviour. I.e., you have to manually turn off agents, not
	// manually turn them on!
	for (bVector::iterator i = eligible.begin(); i != eligible.end(); i++) *i=true;

	agent=0; // Default to first obs, always. Might not be used.

	steps=0;

    }
}
