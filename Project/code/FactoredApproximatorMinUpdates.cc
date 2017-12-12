#include"PGBasics.hh"
#include"BasicController.hh"
#include"FactoredApproximator.hh"
#include"FactoredApproximatorMinUpdates.hh"

namespace libpg {

FactoredApproximatorMinUpdates::FactoredApproximatorMinUpdates(Approximators approximators, bool splitObs)
    : FactoredApproximator::FactoredApproximator(approximators, splitObs) {

    lastUpdateTime.resize(approximators.size());
    hasBeenEligible.resize(approximators.size());
    hasBeenEligibleLongTerm.resize(approximators.size());

    currentTime = 0;
    lastUpdateTime.clear();
    hasBeenEligible.clear();
    hasBeenEligibleLongTerm.clear();
}

void FactoredApproximatorMinUpdates::doApprox(Observation& obs, Vector& output) {
    FactoredApproximator::doApprox(obs, output);
    eligible = &(obs.getEligible());

    for (uint i = 0; i < approximators.size(); i++)
	if ((*eligible)[i]) {
	    hasBeenEligible[i] = true;
	    hasBeenEligibleLongTerm[i] = true;
	}
}


/**
 *
 * The remaining methods just pass onto equiv Approximator
 * functions. We have to do it in this ugly way because more advanced
 * multiple agent controllers will need to override these functions to
 * update all agents in turn.
 *
 */
void FactoredApproximatorMinUpdates::discountTrace() {

    currentTime++;

    int j=0;
    Approximators::iterator i;
    for (i = approximators.begin(); i != approximators.end(); i++) {
	if ((*eligible)(j)) {
	    if (hasBeenEligible(j)) {
		(*i)->setDiscount(pow(discount, currentTime - lastUpdateTime(j)));
		(*i)->discountTrace();
	    }
	    lastUpdateTime(j) = currentTime;
	}
	j++;
    }
}


void FactoredApproximatorMinUpdates::setDiscount(double discount) {
    this->discount = discount;
}


void FactoredApproximatorMinUpdates::instantStep(double reward) {

    for (uint i = 0; i < approximators.size(); i++) {
	if (reward != 0.0) {
	    if (hasBeenEligible(i)) {
		if (!(*eligible)(i)) {
		    // If action eligible, trace has already been updated
		    approximators[i]->setDiscount(pow(discount, currentTime - lastUpdateTime(i)));
		    approximators[i]->discountTrace();
		    lastUpdateTime(i) = currentTime;
		}
		approximators[i]->instantStep(reward);
	    }
	}
    }
}


void FactoredApproximatorMinUpdates::resetTrace() {
    FactoredApproximator::resetTrace();

    currentTime = 0;
    lastUpdateTime.clear();

    /*int k=0;
    for (uint i = 0; i < approximators.size(); i++)
	if (hasBeenEligible(i)) k++;
    std::cerr<<"#eligible actions="<<k;

    k=0;
    for (uint i = 0; i < approximators.size(); i++)
	if (hasBeenEligibleLongTerm(i)) k++;
    std::cerr<<" ("<<k<<")"<<std::endl;
    */
    
    hasBeenEligible.clear();
}

}
