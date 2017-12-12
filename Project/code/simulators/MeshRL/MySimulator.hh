// This is your first RL application. Fill in the blanks shown below
// and in MySimulator.cc. This file is the definition of your control
// problem for reinforcement learning to optimise. You'll probably
// also have to edit TemplateRLApp.cc to get the best performance, but
// your first goal is to fill in the blanks below, compile, and run.

// NOTICE: This template assumes a single RL agent. You should be able
// to figure out what to modify here and in TemplateRLApp.cc and MySimulator.cc
// to make it multi-agent.

// Look for the text [WRITE CODE HERE] to see where to fill stuff in.


#ifndef MySimulator_hh
#define MySimulator_hh

#include <vector>

class MySimulator : public libpg::Simulator {

protected:

    // Insert all the variables and data structures you need
    // to reprsent the state of your problem.
    // [WRITE CODE HERE]

    // Write in your own helper methods if needed.
    // [WRITE CODE HERE]

public:

    // See MySimulator.cc for detailed comments. You don't need to change anything below.

    // Constructor
    MySimulator(char*NodeFileName, char*EleFileName);

    // Empty desctructor. Fill in if you need to deallocate stuff.
    virtual ~MySimulator() {};

    // Not part of the normal simulator interface. Just for TemplateRLApp
    double getDiscountFactor(); 

    virtual void getReward(libpg::Vector& rewards);
    virtual void getObservation(libpg::Observation& obs);
    virtual int doAction(libpg::Vector& action);
    virtual int getMaxEpisodeLength();
    virtual int getObsRows();
    virtual int getObsCols();
    virtual int getAgents();
    virtual int getActionDim();
    virtual int getNumActions();
    virtual int getRewardDim();		
private:
	std::vector<int> nextObtuseVertex;
	int myReward;
	char*mNodeFileName, *mEleFileName;
	size_t num_attempt; // 

};

#endif
