#include"PGBasics.hh"
#include"FactoredApproximator.hh"


using namespace std;

namespace libpg {

    /**
     * Construct a collection of approximators. Implements the many to one controller concept.
     * @param  approximators vector of Approximator* 
     * @param  if true, expect as many colums in the observations as there are approximators
     * and each controller only gets its col as an observation.
     * If false, all approximators get the same observation.
     */
    FactoredApproximator::FactoredApproximator(Approximators approximators, bool splitObs) {

	this->approximators = approximators;
	this->splitObs = splitObs;
	dummyDist.resize(1);

    }


    /**
     * Where the aggregation is done
     */
    void FactoredApproximator::doApprox(Observation& obs, Vector& output) {
    
	// Can only have one winning action here!
	assert(output.size() == approximators.size());

	if (splitObs) {
	    // If we have a col for each approximator, they get their own
	    // observation.
	    assert(obs.getFeatures().size2() == approximators.size());
	} else {
	    assert(obs.getFeatures().size2() == 1);
	}
    
	bool oneEligible = false;

	//cout<<"<<<<<<<<<Approximators.size="<<approximators.size()<<endl;
	// Loop over possible actions, each is its own approximator

	for (unsigned int i = 0; i < approximators.size(); i++) {
	    if (!obs.getEligible()(i)) { // This action not eligible.
		
		// [daa] output[i] = 0.; This is wrong because the number will be exponentiated
		// to make a probability distribution. Doesn't make any difference if the
		// sampler knows about eligible actions, but better not to assume it does.

		output[i] = numeric_limits<double>::min();
	    } else {
		if (splitObs) obs.setAgent(i);
		approximators[i]->doApprox(obs, dummyDist);
		output[i] = dummyDist[0];
		if (isnan(dummyDist[0])) {
		    cout<<"approx "<<i<<" output is nan, weights ";
		    approximators[i]->write(cout);
		    cout<<endl;
		}
		oneEligible = true;
	    }
	}

	if (!oneEligible) throw std::runtime_error("No actions eligible in FactoredApproximator\n");
    }
  

    void FactoredApproximator::feedbackGrad(Observation& obs, Vector& deltas) {

     
	assert(deltas.size() == approximators.size());

	for (uint i=0; i <  approximators.size(); i++) {
	    if (!obs.getEligible()(i)) continue;
	    assert(approximators[i]->getOutputDim() == 1);
	    dummyDist[0] = deltas[i];
	    approximators[i]->feedbackGrad(obs, dummyDist);
	}
	
    }

    
    /**
     *
     * The remaining methods just pass onto equiv Approximator
     * functions. We have to do it in this ugly way because more advanced
     * multiple agent approximators will need to override these functions to
     * update all agents in turn.
     *
     */

    void FactoredApproximator::discountTrace() {

	Approximators::iterator i;
	for (i = approximators.begin(); i != approximators.end(); i++) {
	    (*i)->discountTrace();
	}

    }


    void FactoredApproximator::setDiscount(double discount) {
	Approximators::iterator i;
	for (i = approximators.begin(); i != approximators.end(); i++) {
	    (*i)->setDiscount(discount);
	}
    }


    void FactoredApproximator::accumulateGrad(double reward, Observation& newObs) {

	for (uint i = 0; i < approximators.size(); i++) {
	    if (reward != 0.0) approximators[i]->accumulateGrad(reward, newObs);	
	}
    }


    void FactoredApproximator::resetGrad() {

	Approximators::iterator i;
	for (i = approximators.begin(); i != approximators.end(); i++) {
	    (*i)->resetGrad();
	}

    }


    void FactoredApproximator::instantStep(double reward) {

	for (uint i = 0; i < approximators.size(); i++) {
	    if (reward != 0.0) approximators[i]->instantStep(reward);
	}
    }


    void FactoredApproximator::computeDirection(int steps) { 
	Approximators::iterator i;
	for (i = approximators.begin(); i != approximators.end(); i++) {
	    (*i)->computeDirection(steps);
	}
    }


    void FactoredApproximator::batchStep() { 
	Approximators::iterator i;
	for (i = approximators.begin(); i != approximators.end(); i++) {
	    (*i)->batchStep();
	}
    }

    void FactoredApproximator::setStepSize(double stepSize) { 
	Approximators::iterator i;
	for (i = approximators.begin(); i != approximators.end(); i++) {
	    (*i)->setStepSize(stepSize);
	}
    }


    void FactoredApproximator::resetTrace() { 
	Approximators::iterator i;
	for (i = approximators.begin(); i != approximators.end(); i++) {
	    (*i)->resetTrace();
	}
    }


    void FactoredApproximator::resetParams() { 
	Approximators::iterator i;
	for (i = approximators.begin(); i != approximators.end(); i++) {
	    (*i)->resetParams();
	}
    }


    void FactoredApproximator::randomizeParams(double maxRand) { 
	Approximators::iterator i;
	for (i = approximators.begin(); i != approximators.end(); i++) {
	    (*i)->randomizeParams(maxRand);
	}
    }


    double FactoredApproximator::getMaxParam() {
	double maxi = 0.0;
	Approximators::iterator i;
	for (i = approximators.begin(); i != approximators.end(); i++) {
	    maxi = std::max(maxi, (*i)->getMaxParam());
	}
	return maxi;
    }


    int FactoredApproximator::getNumParams() {
	int p = 0;
	Approximators::iterator i;
	for (i = approximators.begin(); i != approximators.end(); i++) {
	    p += (*i)->getNumParams();
	}
	return p;
    }


    void FactoredApproximator::write(ostream& o) {
	Approximators::iterator i;
	for (i = approximators.begin(); i != approximators.end(); i++) {
	    (*i)->write(o);
	}
    }
  

    void FactoredApproximator::read(istream& o) {
	Approximators::iterator i;
	for (i = approximators.begin(); i != approximators.end(); i++) {
	    (*i)->read(o);
	}
    }  

    int FactoredApproximator::getInputDim() { return approximators[0]->getInputDim(); }
    int FactoredApproximator::getOutputDim() { return approximators.size(); }

    /**
     * @param s what to set, can't do GRADS here.
     */
    void FactoredApproximator::scatter(Vector& v, Approximator::StatsEnum s) {

	int p;
	Approximators::iterator i;

	switch (s) {
        case Approximator::PARAMS:
        case Approximator::GRADS:
        case Approximator::TRACES:
	    p = 0;
	    for (i = approximators.begin(); i != approximators.end(); i++) {
		int nParams = (*i)->getNumParams();
		Vector w(nParams);
		//w.clear();
		w = project(v, range(p ,p+nParams));
		(*i)->scatter(w, s);
		p += nParams;
	    }
	    break;
	default:
	    throw runtime_error("unknown StatsEnum FactoredApproximator::scatter\n");
	}

    }


    void FactoredApproximator::reduce(Vector& v, Approximator::StatsEnum s) {

	int p;
	Approximators::iterator i;

	switch (s) {


	case Approximator::PARAMS:
	case Approximator::GRADS:
	case Approximator::TRACES:
	    p = 0;
	    for (i = approximators.begin(); i != approximators.end(); i++) {
		int nParams = (*i)->getNumParams();
		Vector w(nParams);
		//w.clear();
		(*i)->reduce(w, s);
		project(v, range(p, p+nParams)) = w;
		p += nParams;
	    }

	    break;
	
	default:
	    throw runtime_error("unknown StatsEnum FactoredApproximator::reduce\n");
	}

    }


}
