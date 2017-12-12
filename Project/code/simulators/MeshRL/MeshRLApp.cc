// This is your first RL application. You shouldn't need to edit this
// file until you have MySimulator.cc and MySimlator.hh complete and
// the code compiles and runs. Then you should consider modifying this
// file for best performance and using alternative algorithms. 

// The first thing to consider changing is the algorithm used. This
// file defaults to Natural Actor-Critic. You can also switch to
// Least Squares Policy Iteration. These are both state of the art
// algorithms in RL. See documentation for other altgorithms and
// advanced things to do with TransformControllers.

// The next thing you can do is twiddle the values of the #defines
// for better performance. Note that this code also shows how to
// save parameter files so that you can recover your optimised
// policy.

#include"PGBasics.hh"

// Your simulator
#include"MySimulator.hh"

// Approximators
#include"NeuralNetBatch.hh"
#include"LookupTableBatch.hh" // Use for exact MDPs

// Policies
#include"SoftmaxPolicy.hh"
#include"eGreedyPolicy.hh"

// Controllers
#include"BasicController.hh"
#include"NACTransform.hh"
#include"eGreedyPolicy.hh"
#include"LSTDQController.hh"
#include"ValueController.hh"
#include"SARSAController.hh"
#include"QLearningController.hh"
#include"FactoredController.hh"
#include"LSTDQController.hh"

// RLAlg
#include"OLPomdp.hh"
#include "GOALPomdp.hh"
// Meshing problem implementation 



/**********************************************************
 * Hard coded settings you can change
 **********************************************************/

#define NAC 0
#define LSPI 1
#define ALGORITHM NAC // 0 for Natural actor critic, 1 for LSPI
                      // If you change this, change LAMBDA below.

#define MAX_STEPS 1000000 // Total max optimisation steps, 0=unbounded
#define MAX_TIME 0 // Total max wall clock time, 0=unbounded
#define EPOCH_STEPS 1000 // Steps per epoch (progress print outs/LSTDQ steps)
                         // Increase this if youre getting too much information.
                         // Decrease this if things are very slow, but not
                         // too small for LSPI otherwise it will
                         // increase variance.
                       

// Location of parameter save file
#define PARAMETER_FILE "parameters.txt"

#define LAMBDA 0.9 // For NAC, start with 0.8. Always stay less than
                   // 1.0 Increase slowly only if you're not getting a
                   // good solution. If you're getting no improvement
                   // at all start by increasing the step size below.
                   // For LSPI, start with 0.0 and try values up to an
                   // including 1.0
// LSPI Specific settings
#define EXPLORATION_PROB 0.6 // Probability of a random action

// Natural actor-critic specific settings
#define STEP_SIZE 0.01 // Start small and grow by 10x until things
                          // become unstable during optimisation. If
                          // nothing at all seems to change, this value
                          // might need to be much bigger.
                          // If things improve quickly at the start, then
                          // the reward suddenly gets much worse, make
                          // this smaller.
using namespace std; 
using namespace libpg;


/**
* Temp is a subclass of Temperature class, which
* belongs to SoftmaxPolicy domain. Temp implements
* the user specific temperature decain function.
*/
class Temp : public SoftmaxPolicy::Temperature {

public:

	/**
	* Setting K to 0 makes a constant 1.0 temperature. K>0 implements a temp that
	* decays with inverse square of (steps*K + 1).
	*/
	Temp(double K) {
		this->K = K;
	}

	/**
	* Implements the temperature decaying function.
	* @param steps number of steps done so far.
	* @return temperature value.
	*/
	double getValue(int steps) {
		return (double) 1.0 / ((1.0 + K * (double)steps)*(1.0 + K * (double)steps));
	}

private:

	double K;

};
class EpsilonDecay : public eGreedyPolicy::EpsilonFunction {

private:
	double K;

public:

	EpsilonDecay(double K) {
		this->K = K;
	}

	double getValue(int steps) {
		return (double) 1.0 / ((1.0 + K * (double)steps)*(1.0 + K * (double)steps));
	}

};


RLAlg* createRLAlg_NAC(char*NodeFileName, char*EleFileName) {
    
	/**
	* Natural Actor-Critic controller with linear function approximator	
	**/

    // STEP 1. Create simulator    
	MySimulator* mySimulator = new MySimulator(NodeFileName, EleFileName);
	

    // STEP 2. Create approximator
    Approximator* neuralNet = new NeuralNetBatch(mySimulator->getObsRows(), mySimulator->getNumActions());
    
    // STEP 3. Create a controller 
    Controller* controller;
    
	// Create Natural Actor-Critic controller
	controller = new NACTransform(new BasicController(neuralNet), mySimulator->getDiscountFactor());
      
    // STEP 4. Link all together with an algorithm    
	RLAlg* alg = new GOALPomdp(controller, mySimulator, LAMBDA, STEP_SIZE);
	
    return alg;

}

RLAlg* createRLAlg_SARSA_NeuralNet_Softmax_Decaying(char*NodeFileName, char*EleFileName) {
	
	/**
	* SARSA controller with linear function approximator and
	* Softmax policy (decaying temperature).	
	**/
	double kapa = 1e-05;
	
	// STEP 1. Create simulator 
	MySimulator* mySimulator = new MySimulator(NodeFileName, EleFileName);
	
	// STEP 2. Create approximator
	Approximator* neuralNet = new NeuralNetBatch(mySimulator->getObsRows(), mySimulator->getNumActions());

	// STEP 3. Create a controller 
	Controller* controller = new SARSAController(neuralNet, new SoftmaxPolicy(new Temp(kapa)), mySimulator->getDiscountFactor());

	// STEP 4. Link all together with an algorithm    
	RLAlg* alg = new GOALPomdp(controller, mySimulator, LAMBDA, STEP_SIZE);
	return alg;
}

RLAlg* createRLAlg_SARSA_NeuralNet_eGreedy_Constant(char*NodeFileName, char*EleFileName) {
	
	/**
	* SARSA controller with linear function approximator and
	* e-Greedy policy (constant epsilon).	
	**/
	double epsilon = 0.5;
	// STEP 1. Create simulator 
	MySimulator* mySimulator = new MySimulator(NodeFileName, EleFileName);


	// STEP 2. Create approximator
	Approximator* neuralNet = new NeuralNetBatch(mySimulator->getObsRows(), mySimulator->getNumActions());


	// STEP 3. Create a controller
	Controller* controller = new SARSAController(neuralNet, (Policy*) new eGreedyPolicy(epsilon), mySimulator->getDiscountFactor());

	// STEP 4. Link all together with an algorithm    
	RLAlg* alg = new GOALPomdp(controller, mySimulator, LAMBDA, STEP_SIZE);

	return alg;
}

RLAlg* createRLAlg_SARSA_NeuralNet_eGreedy_Decaying(char*NodeFileName, char*EleFileName) {
	
	/**
	* SARSA controller with linear function approximator and
	* e-Greedy policy (decaying epsilon).	
	**/

	double kapa = 1e-05;
	// STEP 1. Create simulator 
	MySimulator* mySimulator = new MySimulator(NodeFileName, EleFileName);
	
	// STEP 2. Create approximator
	Approximator* neuralNet = new NeuralNetBatch(mySimulator->getObsRows(), mySimulator->getNumActions());

	// STEP 3. Create a controller
	Controller* controller = new SARSAController(neuralNet, (Policy*) new eGreedyPolicy(new EpsilonDecay(kapa)), mySimulator->getDiscountFactor());

	// STEP 4. Link all together with an algorithm    
	RLAlg* alg = new GOALPomdp(controller, mySimulator, LAMBDA, STEP_SIZE);
	return alg;
}

RLAlg* createRLAlg_QLearning_NeuralNet_Softmax_Decaying(char*NodeFileName, char*EleFileName) {
		
	/**
	* QLearning controller with linear function approximator and
	* Softmax policy (decaying temperature).	
	**/

	double kapa = 1e-05;
	// STEP 1. Create simulator 
	MySimulator* mySimulator = new MySimulator(NodeFileName, EleFileName);
	
	// STEP 2. Create approximator
	Approximator* neuralNet = new NeuralNetBatch(mySimulator->getObsRows(), mySimulator->getNumActions());

	// STEP 3. Create a controller
	Controller* controller = new QLearningController(neuralNet, (Policy*) new SoftmaxPolicy(new Temp(kapa)), mySimulator->getDiscountFactor());
	
	// STEP 4. Link all together with an algorithm    
	RLAlg* alg = new GOALPomdp(controller, mySimulator, LAMBDA, STEP_SIZE);
	return alg;
}

RLAlg* createRLAlg_QLearning_NeuralNet_eGreedy_Constant(char*NodeFileName, char*EleFileName) {
	
	/**
	* QLearning controller with linear function approximator and
	* e-Greedy policy (constant epsilon).
	**/

	double epsilon = 0.9;

	// STEP 1. Create simulator 
	MySimulator* mySimulator = new MySimulator(NodeFileName, EleFileName);
		
	// STEP 2. Create approximator
	Approximator* neuralNet = new NeuralNetBatch(mySimulator->getObsRows(), mySimulator->getNumActions());

	// STEP 3. Create a controller
	Controller* controller = new QLearningController(neuralNet, (Policy*) new eGreedyPolicy(epsilon), mySimulator->getDiscountFactor());

	// STEP 4. Link all together with an algorithm    
	RLAlg* alg = new GOALPomdp(controller, mySimulator, LAMBDA, STEP_SIZE);
	return alg;	
}

RLAlg* createRLAlg_QLearning_NeuralNet_eGreedy_Decaying(char*NodeFileName, char*EleFileName) {

	/**
	* QLearning controller with linear function approximator and
	* e-Greedy policy (decaying epsilon).
	**/

	double kapa = 1e-05;

	// STEP 1. Create simulator 
	MySimulator* mySimulator = new MySimulator(NodeFileName, EleFileName);

	// STEP 2. Create approximator
	Approximator* neuralNet = new NeuralNetBatch(mySimulator->getObsRows(), mySimulator->getNumActions());

	// STEP 3. Create a controller
	Controller* controller = new QLearningController(neuralNet, (Policy*) new eGreedyPolicy(new EpsilonDecay(kapa)), mySimulator->getDiscountFactor());

	// STEP 4. Link all together with an algorithm    
	RLAlg* alg = new GOALPomdp(controller, mySimulator, LAMBDA, STEP_SIZE);
	return alg;
}

	

void optimise(char*NodeFileName, char*EleFileName, int algo_num) {
	
	RLAlg* alg = NULL;
	if (algo_num == 1){
		//NAC
		alg = createRLAlg_NAC(NodeFileName, EleFileName);
	}
	else if (algo_num == 2){
		//SARSA with NeuralNet approximator and Softmax policy with decaying temperature
		alg = createRLAlg_SARSA_NeuralNet_Softmax_Decaying(NodeFileName, EleFileName);
	}
	else if (algo_num == 3){
		//SARSA with NeuralNet approximator and e-Greedy policy with constant epsilon
		alg = createRLAlg_SARSA_NeuralNet_eGreedy_Constant(NodeFileName, EleFileName);
	}
	else if (algo_num == 4){
		//SARSA with NeuralNet approximator and e - Greedy policy with decaying epsilon
		alg = createRLAlg_SARSA_NeuralNet_eGreedy_Decaying(NodeFileName, EleFileName);
	}
	else if (algo_num == 5){
		//QLearning with NeuralNet approximator and Softmax policy with decaying temperature
		alg = createRLAlg_QLearning_NeuralNet_Softmax_Decaying(NodeFileName, EleFileName);
	}
	else if (algo_num == 6){
		//QLearning with NeuralNet approximator and e-Greedy policy with constant epsilon
		alg = createRLAlg_QLearning_NeuralNet_eGreedy_Constant(NodeFileName, EleFileName);
	}
	else if (algo_num == 7){
		//QLearning with NeuralNet approximator and e - Greedy policy with decaying temperature
		alg = createRLAlg_QLearning_NeuralNet_eGreedy_Decaying(NodeFileName, EleFileName);
	}
	
    // Indicate to automatically save the parameters whenever performance improves
    // alg->saveBest(PARAMETER_FILE);

    // STEP 5. Learn!    
	alg->learn(EPOCH_STEPS, MAX_TIME, MAX_STEPS);
}

/**
 * Main routine to switch between evaluation and optimisation
 */ 
int main(int argc, char** argv) {

	if (argc == 4){
		if (atoi(argv[3]) == 1){
			printf("\n	Running Actor Critic with NeuralNet approximator");
		}
		else if (atoi(argv[3]) == 2){
			printf("\n	Running SARSA with NeuralNet approximator and Softmax policy with decaying temperature");
		}
		else if (atoi(argv[3]) == 3){
			printf("\n	Running SARSA with NeuralNet approximator and e-Greedy policy with constant epsilon");
		}
		else if (atoi(argv[3]) == 4){
			printf("\n	Running SARSA with NeuralNet approximator and e-Greedy policy with decaying epsilon");
		}
		else if (atoi(argv[3]) == 5){
			printf("\n	Running QLearning with NeuralNet approximator and Softmax policy with decaying temperature");
		}
		else if (atoi(argv[3]) == 6){
			printf("\n	Running QLearning with NeuralNet approximator and e-Greedy policy with constant epsilon");
		}
		else if (atoi(argv[3]) == 7){
			printf("\n	Running QLearning with NeuralNet approximator and e-Greedy policy with decaying epsilon\n");
		}

		for (int i = 0; i < 100; i++){			
			std::cout << "TEST# " << i << std::endl;
			optimise(argv[1], argv[2], atoi(argv[3]));
			std::cout << "**********************************************" << std::endl;
		}
	}
	else{
		printf("\n Invalid Input.");
		printf("\n Usage: $ ./MeshRLApp NodeFileName EleFileName [ALGO_NUM]");
		printf("\n	[ALGO_NUM] could be");
		printf("\n	1 for Actor Critic with NeuralNet approximator");
		printf("\n	2 for SARSA with NeuralNet approximator and Softmax policy with decaying temperature");
		printf("\n	3 for SARSA with NeuralNet approximator and e-Greedy policy with constant epsilon");
		printf("\n	4 for SARSA with NeuralNet approximator and e-Greedy policy with decaying epsilon");
		printf("\n	5 for QLearning with NeuralNet approximator and Softmax policy with decaying temperature");
		printf("\n	6 for QLearning with NeuralNet approximator and e-Greedy policy with constant epsilon");
		printf("\n	7 for QLearning with NeuralNet approximator and e-Greedy policy with decaying temperature\n");

		exit(EXIT_FAILURE);
	}

    exit(EXIT_SUCCESS);

}


