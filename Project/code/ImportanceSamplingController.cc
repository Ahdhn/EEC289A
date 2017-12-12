
// $Id: ImportanceSamplingController.cc 29 2007-01-17 13:17:00Z buffet $

#include"PGBasics.hh"

#include"Sampler.hh"
#include"Teacher.hh"
#include"ImportanceSamplingController.hh"

using namespace std;

namespace libpg {

/**
 * @param  underlying controller.
 */
ImportanceSamplingController::ImportanceSamplingController(Controller* controller, Teacher* teacher, int outputDim)
    : TransformController(controller) {
    this->teacher = teacher;
    isTeaching = true;
    
    distribution.resize(outputDim);
    teacherDistribution.resize(outputDim);
    teacherController = new BasicController(teacher);

    K = 0.0;
    kSamp = 0;
    Ka = 0;
    whichIS = noIS;

}
    

ImportanceSamplingController::~ImportanceSamplingController() {
    delete teacherController;
}


/**
 * Just pass to underlying controller
 */
void ImportanceSamplingController::getAction(Observation& obs, Vector& action, bool computeGrad) {

    if (! isTeaching) {
	// If not observing the teacher, just call the embedded
	// controller's getAction() method.
	
	
	controller->getAction(obs, action, computeGrad);
	
	
    } else {
	/* If observing the teacher: go through the controller's and
	 * the teacher's getAction() methods, so as to compute their
	 * probability distributions over actions.
	 * (plus the gradient for the controller)
	 */
	
	/* 1- Get the distribution produced by the controller,
	 * test if there is an eligible action,
	 * and clean it if a NaN appears.
	 */
	
	// [daa] don't compute the gradient yet, it's the wrong action
	
	controller->getAction(obs, action, false);
	
	distribution.clear();
	controller->reduce(distribution, Approximator::DIST);
	
	Vector distController = distribution;
	
	
	//cout.flush(); cerr.flush();
	
	/* 2- Compute the teacher's distribution. Don't attempt to compute gradients
	 */
	teacherController->getAction(obs, action, false);
	teacherController->reduce(teacherDistribution, Approximator::DIST);
	/* 3- Use these two distributions to compute the distribution to sample from.
	 */
	if (action[0] != -1) {
	    distribution *= (1. - tRate);
	    distribution += tRate*teacherDistribution;
	} else {
	    // Distribution is unchanged.. equals the distController, so K=1
	}
	/*cerr<<"dist(controller+teacher)="<<distribution<<endl;*/
	assert(fabs(norm_1(distribution) - 1.0) < 0.00001);

	/* 4- Pick an action
	 */
	
	int sampledAct  = Sampler::discrete(distribution, 1.0);
	
	action[0] = sampledAct;
	
	
	//cout<<"Action picked="<<action[0]<<endl;
	
	// Update the importance coefficient.
	if ( distController[sampledAct]==0 ||
	     !isnormal(distController[sampledAct]) // Is it a normal floating point number?
	     ) {
	    cout<<"K="<<K<<endl;
	    cout<<"Action picked="<<action[0]<<endl;
	    cout<<distController[sampledAct]<<"/"<<distribution[sampledAct]<<"/";
	    cout<<obs.getEligible()[sampledAct]<<endl;
	    cout<<"dist(controller)="<<distController<<endl;
	    cout<<"dist(controller+teacher)="<<distribution<<endl;
	    exit(0);
	}
	
	double k = 0; 
	
	/* 5- Do the importance sampling re-weighting of the eligibility trace.
	 */
	switch (whichIS) {
	case noIS:
	    break;
	case weightedLocalIS:
	    // [daa] Now that we have settled on an action, compute
	    // the gradient and feedback
	    k = distController[sampledAct] / distribution[sampledAct]; 
	    K += k;
	    kSamp++;
	    break;
	case localOnTraceIS:
	    // Compute importance coefficient, and apply it to the trace.
	    // Can be larger than 1 if distController action chosen is more likey 
	    // according to the FPG controller than the combined distribution. A rare event.
	    // If we get a big reward after this.. then yes, the FPG action should be boosted.
	    K = distController[sampledAct] / distribution[sampledAct];
	    //K = sqrt(distController[sampledAct] / distribution[sampledAct]);
	    controller->setDiscount(K);
	    controller->discountTrace();
	    controller->setDiscount(discount);
	    break;
	default:
	    cout<<"Importance Sampling algorithm not defined !"<<endl;
	    exit(0);
	}
	
	if (computeGrad) controller->computePartialsAndFeedback(obs, sampledAct);
	
	Ka = (1 - 1e-2)*Ka + 1e-2*k;
    }
    
    //cout.flush();
}
    
 /**
  * Just pass to underlying controller
  * @param  Observation vector
  * @param  empty vector to load with the most probably action
  */
 void ImportanceSamplingController::getMostProbAction(Observation& obs, Vector& action) { 

     // Very similar to above method...
     if (! isTeaching)
	 controller->getMostProbAction(obs, action);
     else {
	 controller->getAction(obs, action, false);
	 teacherController->getAction(obs, action);
	 action.clear();

	controller->reduce(distribution, Approximator::DIST);
	teacherController->reduce(teacherDistribution, Approximator::DIST);
	distribution *= (1. - tRate);
	distribution += tRate*(teacherDistribution);

	// Pick most probable action.
	action[0] = 0;
	for (size_t i = 1; i < distribution.size(); i++) {
	    if (distribution[i] > distribution[(int)action[0]]) action[0] = i;
	}
       
    }

}


void ImportanceSamplingController::instantStep(Vector& rewards) {
    /*if (rewards[0]==0 && K==0)
	cout<<"rewards and K are null !\n";*/

    Vector v(rewards.size()); 

    if (!isTeaching)
	controller->instantStep(rewards);
    else {
	switch (whichIS) {
	case noIS:
	case localOnTraceIS:
	    controller->instantStep(rewards);
	    break;
	case weightedLocalIS:
	    K /= kSamp; // Compute averate
	    v = rewards / K; // Divide by average imporance sampling weight
	    controller->instantStep(v);
	    kSamp = 0;
	    K = 0.0;
	    break;
	default:
	    cout<<"Importance Sampling algorithm not defined !"<<endl;
	    exit(0);
	}
    }
}

/*
void ImportanceSamplingController::batchStep() { 
    controller->batchStep(); 
}
*/



void ImportanceSamplingController::setDiscount(double discount) {
    this->discount = discount;
    controller->setDiscount(discount);
}

void ImportanceSamplingController::resetTrace() { 
    //cout<<K<<endl;
    K = 0.0;
    kSamp = 0;
    controller->resetTrace(); 
}


Vector* ImportanceSamplingController::getDistribution() { 
    return (&distribution);
}


void ImportanceSamplingController::turnOffTeaching() {
    if (isTeaching) {
	isTeaching = false;
	delete teacher;
    }
}

void ImportanceSamplingController::setTRate(double d) {
     tRate = d;
}

double ImportanceSamplingController::getTRate() {
    return tRate;
}

void ImportanceSamplingController::setKa(double d) {
    Ka = d;
}

double ImportanceSamplingController::getKa() {
    return Ka;
}

}
