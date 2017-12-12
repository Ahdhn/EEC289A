#ifndef Simulator_hh
#define Simulator_hh
// $Id: Simulator.hh 127 2007-09-10 17:20:04Z daa $

namespace libpg {
/**
 * Container for any RL problem we care to define. We simply have to
 * define the dimensionality of  actions, a reward function, a way to get
 * observationsk, and a way to do actions.
 */

class Simulator {

public:

    /**
     * Shut up compiler warning
     */
    virtual ~Simulator() {};

    /**
     * Maximum episode length in steps. For planning mode. Only used
     * by GOALPomdp at the moment.
     */
    virtual int getMaxEpisodeLength() { return 0; };

    /**
     * Get the number of rows in observation matrix. Generally the
     * length of the observation vector.
     */
    virtual int getObsRows() = 0;

    /**
     * Get the number of columns in the observation matrix Typically
     * use multiple columns to provide a different observation to each
     * agent, so the number of columns is the same as getAgents()
     */
    virtual int getObsCols() = 0;

    /**
     * Get the total number of agents. Usually 1.
     */
    virtual int getAgents() = 0;

    /**
     * Get the total action dimensionality. Often, but not
     * necessarily, the same as the number of agents.
     */
    virtual int getActionDim() = 0;

    /** 
     * Get the number of actions per dimension. Used to initialise the
     * eligible vector in the Observation to the right size
     */
    virtual int getNumActions() = 0;

    /**
     * Get the dimensionality of the reward vector. Often, but not
     * necessarily, the same as the number of agents.
     */
    virtual int getRewardDim() = 0;

    /**
     * Just return the reward associated with the most recent state achieved.
     * @param  vector of rewards. The vector has getAgents() dimensionality
     * In non-local reward mode this vector will just be averaged.
     */
    virtual void getReward(Vector& rewards) = 0;

    /**
     * Return the observation of the current state. Could be a scalar,
     * vector, or matrix, so we left the type as Matrix
     * @return  possibly stochastic observation of the current state
     */
    virtual void getObservation(Observation& obs) = 0;

    /**
     * Do action issued by a controller @return  vector of actions,
     * typically just a scalar which is cast to an int, but it there
     * are multiple agents it will be a vector of actions for each
     * agent.
     * @param  vector that describes the action to perform.
     */
    virtual int     doAction(Vector& action) = 0;


};
}
#endif
