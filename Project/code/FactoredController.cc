#include"PGBasics.hh"
#include"FactoredController.hh"

using namespace std;
namespace libpg {

/**
 * @param controllers vector of Controller* 
 * @param splitObs if true,
 * expect as many colums in the observations as there are controllers
 * and each controller only gets its row as an observation. If false,
 * all controllers get the same observation.
 * @param localRewards will each controller end up getting its own reward, i.e., does Factored controller split up the reward vector
 */
FactoredController::FactoredController(Controllers controllers, bool splitObs, bool localRewards) {

    this->controllers = controllers;
    this->splitObs = splitObs;
    this->localRewards = localRewards;
    dummyAction.resize(1);
    subReward.resize(1);
    
}


/**
 * Where the real coding is done to implement softmax
 */
void FactoredController::getAction(Observation& obs, Vector& action,  bool computeGrad) {

    // One action for each controller
    assert(action.size() == controllers.size());

    if (splitObs) {
	// If we have a col for each controller, they get their own
	// observation.
	assert(obs.getFeatures().size2() == controllers.size());
    } else {
	assert(obs.getFeatures().size2() == 1);
    }
	
    action.clear(); // Actions aren't necessarily cleared at higher
		    // level. Not an issue for other styles of
		    // controller
    for (unsigned int i = 0; i < controllers.size(); i++) {
 	if (!obs.getEligible()[i]) continue; // This action not eligible. 	
	if (splitObs) obs.setAgent(i);
	controllers[i]->getAction(obs, dummyAction, computeGrad);
	action[i] = dummyAction[0];
    }
}


/**
 * Just like above, but always returns each controller most probable action
 * @param  observation features
 * @param  empty array of actions to populate
 */
void FactoredController::getMostProbAction(Observation& obs, Vector& action) {

    // One action for each controller
    assert(action.size() == controllers.size());

    if (splitObs) {
	// If we have a col for each controller, they get their own
	// observation.
	assert(obs.getFeatures().size2() == controllers.size());
    } else {
	assert(obs.getFeatures().size2() == 1);
    }
	
    action.clear(); // Actions aren't necessarily cleared at higher
		    // level. Not an issue for other styles of
		    // controller
    for (unsigned int i = 0; i < controllers.size(); i++) {
	if (!obs.getEligible()[i]) continue; // This action not eligible.
	if (splitObs) obs.setAgent(i);
	controllers[i]->getMostProbAction(obs, dummyAction);
	action[i] = dummyAction[0];
    }
}


int FactoredController::getInputDim() {
    return controllers[0]->getInputDim();
}


/**
 *
 * The remaining methods just pass onto equiv Approximator
 * functions. We have to do it in this ugly way because more advanced
 * multiple agent controllers will need to override these functions to
 * update all agents in turn.
 *
 */

void FactoredController::discountTrace() {

    Controllers::iterator i;
    for (i = controllers.begin(); i != controllers.end(); i++) {
	(*i)->discountTrace();
    }

}


void FactoredController::setDiscount(double discount) {
    Controllers::iterator i;
    for (i = controllers.begin(); i != controllers.end(); i++) {
	(*i)->setDiscount(discount);
    }
}


void FactoredController::accumulateGrad(Vector& rewards, Observation& newObs) {


    assert((localRewards && rewards.size()==controllers.size()) || 
	   (!localRewards && rewards.size()==1));
    

    for (uint i = 0; i < controllers.size(); i++) {
	subReward[0] = rewards[localRewards?i:0];
	controllers[i]->accumulateGrad(subReward, newObs);
    }
}


void FactoredController::resetGrad() { 
    Controllers::iterator i;
    for (i = controllers.begin(); i != controllers.end(); i++) {
	(*i)->resetGrad();
    }
}


void FactoredController::instantStep(Vector& rewards) {


   assert((localRewards && rewards.size()==controllers.size()) || 
	  (!localRewards && rewards.size()==1));


    for (uint i = 0; i < controllers.size(); i++) {
	subReward[0] = rewards[localRewards?i:0];
	controllers[i]->instantStep(subReward);
    }
}
/*
void FactoredController::instantStep(Vector& rewards, Observation& newObs) {


   assert((localRewards && rewards.size()==controllers.size()) || 
	  (!localRewards && rewards.size()==1));


    for (uint i = 0; i < controllers.size(); i++) {
	subReward[0] = rewards[localRewards?i:0];
	controllers[i]->instantStep(subReward, newObs);
    }
}
*/

void FactoredController::computeDirection(int steps) { 
    Controllers::iterator i;
    for (i = controllers.begin(); i != controllers.end(); i++) {
	(*i)->computeDirection(steps);
    }
}


void FactoredController::batchStep() { 
    Controllers::iterator i;
    for (i = controllers.begin(); i != controllers.end(); i++) {
	(*i)->batchStep();
    }
}


void FactoredController::setStepSize(double stepSize) { 
    Controllers::iterator i;
    for (i = controllers.begin(); i != controllers.end(); i++) {
	(*i)->setStepSize(stepSize);
    }
}


void FactoredController::resetTrace() { 
    Controllers::iterator i;
    for (i = controllers.begin(); i != controllers.end(); i++) {
	(*i)->resetTrace();
    }
}


void FactoredController::resetParams() { 
    Controllers::iterator i;
    for (i = controllers.begin(); i != controllers.end(); i++) {
	(*i)->resetParams();
    }
}


void FactoredController::randomizeParams(double maxRand) { 
    Controllers::iterator i;
    for (i = controllers.begin(); i != controllers.end(); i++) {
	(*i)->randomizeParams(maxRand);
    }
}


double FactoredController::getMaxParam() {
    double maxi = 0.0;
    Controllers::iterator i;
    for (i = controllers.begin(); i != controllers.end(); i++) {
	maxi = std::max(maxi, (*i)->getMaxParam());
    }
    return maxi;
}


int FactoredController::getNumParams() {
    int p = 0;
    Controllers::iterator i;
    for (i = controllers.begin(); i != controllers.end(); i++) {
	p += (*i)->getNumParams();
    }
    return p;
}


void FactoredController::write(ostream& o) {
    Controllers::iterator i;
    for (i = controllers.begin(); i != controllers.end(); i++) {
	(*i)->write(o);
	o<<endl;
    }
}
  

void FactoredController::read(istream& o) {
    Controllers::iterator i;
    for (i = controllers.begin(); i != controllers.end(); i++) {
	(*i)->read(o);
    }
}  


void FactoredController::reduce(Vector& v, Approximator::StatsEnum s) {
    Controllers::iterator i;
    Vector subV;
    int p = 0;

    assert(v.size() == (size_t)getNumParams());

    for (i = controllers.begin(); i != controllers.end(); i++) {
	subV.resize((*i)->getNumParams());
	(*i)->reduce(subV, s);
	project(v, range(p, p+subV.size())) += subV;
	p += subV.size();
    }
}


void FactoredController::scatter(Vector& v, Approximator::StatsEnum s) {
    Controllers::iterator i;
    Vector subV;
    int p = 0;

    assert(v.size() == (size_t)getNumParams());

    for (i = controllers.begin(); i != controllers.end(); i++) {    
	subV.resize((*i)->getNumParams());
	subV.assign(project(v, range(p, p+subV.size())));
	(*i)->scatter(subV, s);
	p += subV.size();
    }
}
}
