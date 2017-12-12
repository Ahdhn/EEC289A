// This is your first RL application. Fill in the blanks shown below
// and in MySimulator.hh. This file is the definition of your control
// problem for reinforcement learning to optimise. You'll probably
// also have to edit TemplateRLApp.cc to get the best performance, but
// your first goal is to fill in the blanks below, compile, and run.

// NOTICE: This template assumes a single RL agent. You should be able
// to figure out what to modify here and in TemplateRLApp.cc to make
// it multi-agent.

// Look for the text [WRITE CODE HERE] to see where to fill stuff in.

#include"PGBasics.hh"
#include"MySimulator.hh"

#include "src/MeshOptviaSampling2D_mini.hh"

#include <iostream>
#include <vector>
#include <algorithm> 
// You can use normal C++ symbols
using namespace std;

// You can access the libpg library
using namespace libpg;

/**
 * Constructor for your simulator. Do any once off
 * initialisation. Read data files, allocate memory, stuff like that.
 * Can be left empty if need be.
 */
MySimulator::MySimulator(char*NodeFileName, char*EleFileName) {

    // [WRITE CODE HERE]

	srand (time(NULL) ); // activate for different experiments
	num_attempt = 0;
	InitiateRandNumGenerator(rand());
	mNodeFileName = NodeFileName;
	mEleFileName = EleFileName;
	ReadTriangleInput(NodeFileName, EleFileName);
	
	VerifyMeshNonObtuse(_min_angle_input, _num_obtuse_input);
	
	//std::cout << "Input number of samples = " << _num_samples << std::endl;
	//std::cout << "Input Min angle = " << _min_angle_input << std::endl;
	//std::cout << "Input Max angle = " << _max_angle_input << std::endl;
	//std::cout << "Input number of obtuse angles = " << _num_obtuse_input << std::endl;


	//create a vector with random entry from which we gonna pick the next obtuse triangles 
	//searching in it
	nextObtuseVertex.resize(_num_samples);
	for (size_t i = 0; i < _num_samples; i++){
		nextObtuseVertex[i] = i;
	}
	std::random_shuffle(nextObtuseVertex.begin(), nextObtuseVertex.end());

	
	//NonObtuseReTriangulation(_min_angle_input);	
	//double max_ang, min_ang;
	//get_max_min_angles(_num_samples, _del, _samples, max_ang, min_ang, FindApex);
	//std::cout << "\n" << std::endl;
	//size_t num_obtuse;
	//VerifyMeshNonObtuse(_min_angle_input, num_obtuse);	
	//std::cout << "Output number of samples = " << _num_samples << std::endl;
	//std::cout << "Output Min angle = " << min_ang << std::endl;
	//std::cout << "Output Max angle = " << max_ang << std::endl;
	//std::cout << "Output number of obtuse angles = " << num_obtuse << std::endl;
	
	//plot_unit_box("input.ps", _samples, _del, _num_samples, 0.02, 1, 0, 0, 1, 0, 0, 0, 0, NULL, 0, 0, NULL);
}



/**
 * Return the reward for the current state, or previous-state and action combination.
 * In the multi-agent case there is one reward per reward. Here we assume a single agent, so only
 * a single scalar reward. It's often easier to compute the reward during doAction() and save
 * it a class variable to return here.
 * @param rewards vector to put the rewards in. 
 */
void MySimulator::getReward(Vector& rewards) {

    // Check single agent assumption
    assert(rewards.size() == 1);

    // Initialise reward
    //double reward = 0;

	// Insert scalar reward into return vector
    // [WRITE CODE HERE]
    // reward = ...
	
    //rewards[0] = reward;		
	rewards[0] = myReward;

}


/**
 * Fill in a description of the current state. Think of the
 * observation for a single agent as being a single column matrix of
 * dimension getObsRows()*1. 
 *
 * An observation is sometimes refered to as a basis feature.
 *
 * It's up to you to define below in
 * getObsRows() how big the observation is. Making it really big
 * can degrade performance because more parameters are required.
 * Choose your observation features carefully.  
 *
 * One tip is to try and normalise each element of your observation
 * to be between -1 and 1. Another tip is to make sure at least
 * one element of the observation is a constant 1.0 to provide
 * bias to the neural network approximator (don't worry if you don't
 * know what this means, just do it.)
 *
 * In an MDP setting there would typically be only one
 * element in the observation, the state-ID. This is not
 * appropriate for TemplateRLApp. Use a state basis feature
 * instead. If you want to use state-IDs for observations, edit
 * TemplateRLApp.cc to use TableLookup controllers instead of
 * NeuralNet controllers.
 *
 * @param obs The observation class to fill in and return
 */
void MySimulator::getObservation(Observation& obs) {

    // Retrieve the feature matrix from the observation.
    Matrix myObs = obs.getFeatures();

    // Check size assumptions
    assert(myObs.size1() == (size_t)getObsRows());
    assert(myObs.size2() == 1); // Single agent

    myObs.clear();

    // Fill in the observations for each row of the matrix
    // There's only one column, column 0.
    // [WRITE CODE HERE]
    // myObs(0, 0) = some value;
    // myObs(1, 0) = some value;
    // ...
    // myObs(getObsRows()-1, 0) = some value;

    // Note, don't do anything that affects the state of the 
    // simulator. You never know if getObservation() will
    // be called more than once between doAction()
	
	//we observe the following
	//# obtuse triangles / # total triangles 
	//# successful ejection / # failed ejection
	//# successful attractor ejection / # failed attractor ejection
	//# successful injection / # failed injection
	//# successful repeller injection / # failed repeller injection

	myObs(0, 0) = NormalizeNumObtuseTriangles();
	//for (size_t ip = 1; ip <= _num_samples_input; ip++){
	//	if (ip < _num_samples){
	//		myObs(2 * ip - 1, 0) = _samples[ip - 1][0];
	//		myObs(2 * ip, 0) = _samples[ip - 1][1];
	//	}
	//	else{
	//		myObs(2 * ip - 1, 0) = 0;
	//		myObs(2 * ip, 0) = 0;
	//	}
	//}

	//myObs(1, 0) = double(_numSifting) / double(_numSifting_failed);
	//myObs(2, 0) = double(_numSiftingAtt) / double(_numSiftingAtt_failed);
	//myObs(3, 0) = double(_numInject) / double(_numInject_failed);
	//myObs(4, 0) = double(_numInjectRep) / double(_numInjectRep_failed);
}


/**
 * Do the action passed in. Again, in the single agent case the action vector
 * will always be 1x1, and contain an integer valued action.
 * You should take this action and update the state of the simulator 
 * accordingly. You can add state variables in MySimulator.hh
 * IMPORTANT: For this template I assume an infinite-horizon process.
 * This means the simulation does not terminate. If your problem
 * is finite-horizon (it has terminal state or goals), then it's your job
 * to reset the state to the initial state when you encounter a terminal 
 * state here. By editing TemplateRLApp.cc you can choose GoalPOMDP as the
 * top level algorithm appropriate for episodic cases.
 * @param action The 1x1 vector containing the action to do.
 * @return 1 for goal state (ignored in TemplateRLApp).
 */
int MySimulator::doAction(Vector& action) {

    // Check the assumption on the size of the action
    assert(action.size() == 1);
   
    int actionID = (int)action[0]; // Will generate a compiler warning until used
	
    // [WRITE CODE HERE]
    // Do stuff with actionID to update the state

	
	/********************DO THE ACTION********************/
	bool success = false;
	if (actionID == 0){
		//relocate all vertices 
		for (size_t ip = 0; ip<_num_samples; ip++){
			NonObtuseRelocate_Caller(ip, _min_angle_input);
		}
		success = true;
	}
	else if (actionID == 1){
		//ejecte one vertex 
		size_t n1, n2;
		std::random_shuffle(nextObtuseVertex.begin(), nextObtuseVertex.end());
		for (size_t i = 0; i < nextObtuseVertex.size(); i++){
			//check if this one has an obtuse angle and send it 
			//for ejection if it is true and then quit 
			if (!ObtuseHead(nextObtuseVertex[i], n1, n2)){ continue; }
	
			if (NonObtuseSift_Caller(nextObtuseVertex[i], n1, n2, _min_angle_input)){
				_numSifting++;
				success = true;
			}
			else{
				_numSifting_failed++;
			}
			break;
		}		
	}
	else if (actionID == 2){
		//eject attractor one vertex 
		size_t n1, n2;
		std::random_shuffle(nextObtuseVertex.begin(), nextObtuseVertex.end());
		for (size_t i = 0; i < nextObtuseVertex.size(); i++){
			//check if this one has an obtuse angle and send it 
			//for attractor ejection if it is true and then quit 
			if (!ObtuseHead(nextObtuseVertex[i], n1, n2)){ continue; }
			if (NonObtuseSiftAttractor_Caller(nextObtuseVertex[i], n1, n2, _min_angle_input)){
				_numSiftingAtt++;
				success = true;
			}
			else{
				_numSiftingAtt_failed++;
			}
			break;
		}
	
	}
	else if (actionID == 3){
		//inject one vertex 
		size_t n1, n2;
		std::random_shuffle(nextObtuseVertex.begin(), nextObtuseVertex.end());
		for (size_t i = 0; i < nextObtuseVertex.size(); i++){
			//check if this one has an obtuse angle and send it 
			//for injection if it is true and then quit 
			if (!ObtuseHead(nextObtuseVertex[i], n1, n2)){ continue; }
			if (NonObtuseInjection(nextObtuseVertex[i], n1, n2, 0, _min_angle_input)){
				_numInject++;
				success = true;
			}
			else{
				_numInject_failed++;
			}
			break;
		}
	}
	else if (actionID == 4){
		//inject repel one vertex 
		size_t n1, n2;
		std::random_shuffle(nextObtuseVertex.begin(), nextObtuseVertex.end());
		for (size_t i = 0; i < nextObtuseVertex.size(); i++){
			//check if this one has an obtuse angle and send it 
			//for repeller injection if it is true and then quit 
			if (!ObtuseHead(nextObtuseVertex[i], n1, n2)){ continue; }
			if (NonObtuseInjectRepeller(nextObtuseVertex[i], n1, n2, 0, _min_angle_input)){
				_numInjectRep++;
				success = true;
			}
			else{
				_numInjectRep_failed++;
			}
			break;
		}
	
	}

	
	if (success && (actionID == 1 || actionID == 2)){
		//remove _num_samples from nextObtuseVertex
		nextObtuseVertex.erase(std::remove(nextObtuseVertex.begin(), nextObtuseVertex.end(), _num_samples), nextObtuseVertex.end());
	}
	else if (success && (actionID == 3 || actionID == 4)){
		//insert the new added vertex and re-shuffle 
		nextObtuseVertex.push_back(_num_samples - 1);
		std::random_shuffle(nextObtuseVertex.begin(), nextObtuseVertex.end());
	}


	/********************GET REWARD********************/
	//calc the reward now and store it in myReward	
	VerifyMeshNonObtuse(_min_angle_input, _current_num_non_obtuse);
	num_attempt++;
	if (_current_num_non_obtuse == 0){
		//goal reached 
		myReward = 1000;
		ResetMesh(mNodeFileName, mEleFileName);
		num_attempt = 0;
	}
	else if (success){ 
		//if it was successful in eliminating an obtuse angle 
		myReward = 1;
	}
	else if (!success){
		myReward = 0;
	}
	if (_current_num_non_obtuse != 0 && num_attempt == 10 * _num_samples_input){
		//if we reached the max number of attempts without having
		//non-obtuse meshing, give penality and reset mesh
		//because it means we are stucked 
		//std::cout << " STUCK :( " << std::endl;
		myReward = -100;
		
		//plot_unit_box("stuck.ps", _samples, _del, _num_samples, 0.02, 1, 0, 0, 1, 0, 0, 0, 0, NULL, 0, 0, NULL);

		ResetMesh(mNodeFileName, mEleFileName);
		
		_current_num_non_obtuse = 0;

		num_attempt = 0;
	}
	
	//std::cout << "num_attempt= " << num_attempt << " myReward" << myReward << std::endl;
	//return 0;
	return  _current_num_non_obtuse == 0; // This indicates goal state not encountered.
                                          // Always return 0 for infinite horizon.
                                          // This actually won't be checked by TemplateRLApp
                                          // Unless you change which RLAlg is used.
}



/**
 * This defines how many elements there are in your observation vector.
 * This can be as big as you want. Minimum 1 element.
 * 
 * If you in an MDP setting, you can just say getObsRows=1 and you 
 * put the state-id in that element. But see warning above.... you
 * need to edit TempateRLApp.cc to do this properly.
 */
int MySimulator::getObsRows() {
  
    int obsRows;

    // [WRITE CODE HERE] 
    // obsRows = how many elements are there per observation
		
	//obsRows = 5;
	obsRows = 1;	
	//obsRows = 2 * _num_samples_input + 1;
    return obsRows;
}


/**
 * Return the number of actions the simulator knows about.
 * The controller will only produce actionID's between
 * 0 and getNumActions()-1
 */
int MySimulator::getNumActions() {

    int numActions;

    // [WRITE CODE HERE]
    // numActions = how many actions do you have?
	
	//We have five actions
	//1 relocate all vertices (smothing)-> too expensive 
	//2 ejection one vertex -> favourite 
	//3 ejection attractor one vertex -> second favourite 
	//4 injection one vertex -> okay
	//5 inection repel one vertex -> not bad
	numActions = 5; 

    return numActions;
}



/**
 * The discount factor for this problem. If you don't know what a discount factor is
 * then read Barto and Sutton "Reinforcement Learning: An introduction", 1998.
 * Note this is not part of the standard simulator interface, but many RL people
 * expect the discount factor to be defined as part of the problem. For infinite
 * horizon settings (assumed by TemplateRLApp) this needs to be between [0,1).
 * Do not choose 1.0!
 */
double MySimulator::getDiscountFactor() {
		
    double discount = 0.9;

    // [WRITE CODE HERE]
    // discount =

    assert(discount >= 0.0);
    assert(discount < 1.0);
    return discount;

}


/***************************************************************
 * You don't need to edit any lower than this line for the
 * single-agent case.
 ***************************************************************/


/**
 * Get the dimensionality of the reward.
 * Hard wired to 1 for single agent case.
 * Only different if multi-agents and each
 * agent generates its own local reward.
 */
int MySimulator::getRewardDim() {		
    return 1;
}


/**
 * Get the total number of agents. Hard wired to 1 here.
 */
int MySimulator::getAgents() {
    return 1;
}


/** 
 * Maximum episode length in steps. For planning mode. Only used
 * by GOALPomdp at the moment.
 * Do not change.
 */
int MySimulator::getMaxEpisodeLength() {
	return 10 * _num_samples_input;
}


/**
 * Get the number of columns in the observation matrix Typically
 * use multiple columns to provide a different observation to each
 * agent, so the number of columns is the same as getAgents()
 */
int MySimulator::getObsCols() {
    return getAgents();
}


/**
 * Get the total action dimensionality. Often, but not
 *  necessarily, the same as the number of agents.
 */
int MySimulator::getActionDim() {
    return getAgents();
}


