#ifndef Policy_hh
#define Policy_hh

namespace libpg {
/**
 * This class defines the policy used to select the next action.
 * Each policy will be implemented as a subclass of this kind.
 */

class Policy {


public:

    /**
     * Shut up compiler warning.
     */
    virtual ~Policy() {}; 

    /**
     * Sample next action according to policy, given a vector and
     * observations. Note that the dist value vector is passed by copy deliberately.
     * @param dist the vector of values for each action
     * @param action return value indicting appropriate action(s)
     * @param obs the observation vector
     */
    virtual void getAction(Vector dist, Vector& action, Observation& obs) = 0;


};
}
#endif
