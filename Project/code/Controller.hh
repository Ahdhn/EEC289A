#ifndef Controller_hh
#define Controller_hh
// $Id: Controller.hh 127 2007-09-10 17:20:04Z daa $


namespace libpg {
/**
 * A controller is an agent that represents the parameterised stochastic action 
 * selector. Apart from choosing actions, it also stores parameters, 
 * computes gradients, and knows about the basic GPOMDP eligibility trace 
 * update maths and gradient step maths.
 *
 * This is just the abstract version with necessary methods
 */
class Controller {

public:
   
    /**
     * Shut up compiler warning
     */
    virtual ~Controller() {};

    /**
     * Main work horse, feed observations through, sample an action.
     * Optionally accumulate instant log action grad directly to trace
     */
    virtual void getAction(Observation& obs, Vector& action, bool computeGrad = true) = 0;


    /**
     * compute partial derivative of prob of action wrt a vector of likelihoods
     * then feedback to the underlying parametersiation.
     */
    virtual void computePartialsAndFeedback(Observation& obs, int sampledAct, double weight = 1.0) = 0;

    /**
     * Instead of sampling an action, pick the one with the highest prob.
     * Use to exploit policy after learning.
     */
    virtual void getMostProbAction(Observation&, Vector&) { NYI }


    /**
     * Find out what the input dimension of this controller is
     */
    virtual int getInputDim() = 0;

    /**
     * How many parameters does it have
     */
    virtual int getNumParams() = 0;

    /**
     * Remaining methods just pass through to equivalent in 
     * Approximator
     */
    virtual void discountTrace() = 0;

    // The lambda value
    virtual void setDiscount(double) = 0;
    virtual void accumulateGrad(Vector&, Observation&) { NYI };
    virtual void resetGrad() { NYI };
    virtual void computeDirection(int) { NYI };
    virtual void instantStep(Vector&) {NYI};
    virtual void batchStep() { NYI };
    virtual void setStepSize(double) = 0;
    virtual void resetTrace() = 0;
    virtual void resetParams() = 0;

    virtual void randomizeParams(double) { NYI };
    virtual double getMaxParam() = 0;
    virtual void write(std::ostream&) { NYI };
    virtual void read(std::istream&) { NYI };


    /**
     * Add either the parameters or gradient or traces into a vector 
     * @param v vector to store stuff into. v must be at least as big
     * as total params/grads/traces. Will not be resized. Will not be cleared.
     * @param s whether to collect grads or params or traces
     */
    virtual void reduce(Vector&, Approximator::StatsEnum) { NYI };
    
    
    /**
     * Disperse either the parameters or gradient into appropriate
     * data structures for this class.
     * @param v Vector to replace current params or grads or
     * traces. Must be correct size
     * @param  s whether to restore grads or params or traces
     */
    virtual void scatter(Vector&, Approximator::StatsEnum) { NYI };

  
};
}
#endif
