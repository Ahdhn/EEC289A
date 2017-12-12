#ifndef Observation_hh
#define Observation_hh
// $Id: Observation.hh 127 2007-09-10 17:20:04Z daa $

namespace libpg {
class Observation {

protected:



    /**
     * The actual full observation, filled in by
     * the simulator
     */
    Matrix features;
    
    /**
     * A vector that describes which action are available It's up to a
     * controller to decide how to use this, if at all. It is initialised
     * to dimension=number of agents if agents>1 and the number of actions
     * otherwise.
     *
     * Note: This is not a boost::numeric::ublas vector, but an STL one !
     */
    bVector eligible;
    
    /**
     * Another optional variable designed for multiple controllers
     * where the observation for each agent is encoded as a row
     * or column of the observation vector. Determines which
     * row or col should be used in the current controller.
     * Basically an efficiency thing, avoiding copying a 
     * potentially large observation column into a stand alone
     * column
     */
    int agent; 

    /**
     * A count of the number of steps in optimisation/evaluation so
     * far.  This may only be reset in between epochs (in the case of
     * batch learning) Or after a complete learning reset (i.e., reset
     * of parameters).  Any part of the controller/agent code may use
     * this. E.g., if the policy is non-stationary, or if step sizes
     * need to be computed based on optimisation steps.  It may only
     * be set by the RLAlg class.
     */
   int steps;

public:

    Observation(int rows, int cols, int agents, int actions);

    Observation(Observation& obs);
    Observation() {}; // Allow creation of a blank Observation

    /**
     * Resize the observation matricies and vectors.
     * Normally done by the constructor with the same arguments.
     * Elgible vector size inited to agents if agents>1, actions otherwise 
     */
    void init(int rows, int cols, int agents, int actions);
    
    /**
     * Get count of the number of steps in optimisation/evaluation so
     * far.  This may only be reset in between epochs (in the case of
     * batch learning) Or after a complete learning reset (i.e., reset
     * of parameters).  Any part of the controller/agent code may use
     * this. E.g., if the policy is non-stationary, or if step sizes
     * need to be computed based on optimisation steps.  It may only
     * be set by the RLAlg class.
     */
    int getSteps() { return steps; }

    /**
     * Set the number of steps. Only RLAlg should call this.
     * @param s new step number
     */
    void setSteps(int s) { steps = s; }


    /**
     * set which agent this observation is destined for. Used
     * by approximators and controllers that single out a particular
     * part of the feature matrix depending on the agent.
     */
    void setAgent(int a) { agent = a; }

    /**
     * Find out which agent is being represented now. Assists with
     * pulling out the correct part of the features, but exactly
     * how this is done depends on the controller and approximators 
     * involved.
     */
    int getAgent() { return agent; }

    /**
     * Get the observation matrix.
     */
    Matrix& getFeatures() { return features; }

    /**
     * Get the eligibility vector
     */
    bVector& getEligible() { return eligible; }
    
    /**
     * Perfect copy of observation, including resizing
     */
   Observation& operator=(Observation& rhs);

};
}
#endif
