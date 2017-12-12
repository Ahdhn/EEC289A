Installation is done on Linux environment 

Dependencies:    
    Boost 
    uBlas sandbox
    ATLAS
    CLAPACK
	
Installation:
	1- Edit Makefile to set appropriate paths and defines for optional packages. If everything is set to 
	the default location, no edit is needed. 	
	2- '$ make' will give you libpg.a
	3- Go into simulators/MeshRL
	4- Again, edit Makefile to the libraries
    5-'$ make' will give you  working ./MeshRLApp
	6- Use the following command to run the modle inside /input 
	  '$./MeshRLApp ./MeshRLApp input/dh_102_nodes_minA_35.09.node input/dh_102_nodes_minA_35.09.ele 2'
	Use different number (at the end) for different algorithm 
	   1 for Actor Critic with NeuralNet approximator
	   2 for SARSA with NeuralNet approximator and Softmax policy with decaying temperature
	   3 for SARSA with NeuralNet approximator and e-Greedy policy with constant epsilon
	   4 for SARSA with NeuralNet approximator and e-Greedy policy with decaying epsilon
	   5 for QLearning with NeuralNet approximator and Softmax policy with decaying temperature
	   6 for QLearning with NeuralNet approximator and e-Greedy policy with constant epsilon
	   7 for QLearning with NeuralNet approximator and e-Greedy policy with decaying temperature