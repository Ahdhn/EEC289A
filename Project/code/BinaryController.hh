#ifndef BinaryController_hh
#define BinaryController_hh
// $Id: BinaryController.hh 127 2007-09-10 17:20:04Z daa $

namespace libpg {
/**
 * Identical to BasicController except that the output dimension of
 * the approximator must be 1, and we output the action as either 1,
 * or 0 based on a logistic regression mapping of the single approx
 * output to a probability of choosing 1 or 0. Halves the number of
 * parameters needed in the case where decisions are binary.
 */

class BinaryController : public BasicController {

public:
    
    BinaryController(Approximator* approx);

    /**
     * Where the real coding is done to implement 
     */
    virtual void getAction(Observation& obs, Vector& action, bool computeGrad = true);
    
    virtual void getMostProbAction(Observation& obs, Vector& action);

    
};
}
#endif
