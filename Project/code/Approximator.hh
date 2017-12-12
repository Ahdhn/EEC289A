#ifndef Approximator_hh
#define Approximator_hh
// $Id: Approximator.hh 91 2007-04-11 06:05:23Z daa $

namespace libpg {

/**
 * Underlying parameterised function approximator that approximated
 * distributions, either discrete or continuous.
 */

/**
 * Knows about, and performs, basic policy, or value, gradient operations.
 */

    class Approximator {
    
    
    public:
	
	/**
	 * Shut up compiler warning
	 */
	virtual ~Approximator() {};
	
	
	/**
	 * Allows operations that could be equally be performed on params
	 * or grads or traces to choose which.
	 */
	enum StatsEnum { PARAMS, GRADS, TRACES, DIST };
	
	/**
	 * Constructor, do whatever you need to.
	 * Approximator() {};
	 */
	
	virtual int getNumParams() = 0;
	
	/**
	 * Main work horse, feed observations through to get the distribution 
	 * of actions
	 */
	virtual void doApprox(Observation& obs, Vector& output) = 0;
	
	/**
	 * Feed back a gradient signal. Final gradient accumulated directly to
	 * trace object in policy.
	 */
	virtual void feedbackGrad(Observation& obs, Vector& deltas) = 0; 
	
	/**
	 * Just multiply trace by  discount factors
	 */
	virtual void discountTrace() = 0;
	
	/**
	 * Set the discount factors. Does not necessarily have to be constant.
	 */
	virtual void setDiscount(double) = 0;
	
	/**
	 * Calculate instant gradient and step the parameters
	 */
	virtual void instantStep(double) = 0;
	
	/**
	 * Set the step size. Does not necessarily have to be a constant.
	 */
	virtual void setStepSize(double) = 0;
	
	/**
	 * Reset the eligibility trace to 0
	 */
	virtual void resetTrace() = 0;
	
	/**
	 * Do whatever resetting the params means.
	 */
	virtual void resetParams() = 0;
	
	/**
	 * Create a new random set of params.
	 */
	virtual void randomizeParams(double) { NYI };
	
	/**
	 * Retrieve largest absolute parameter value.
	 */
	virtual double getMaxParam() = 0;
	
	/**
	 * Read or write parameters
	 */
	virtual void write(std::ostream&) { NYI };
	virtual void read(std::istream&) { NYI };
	
	
	/**
	 * Query expected number of inputs
	 */
	virtual int getInputDim() = 0;
	
	/**
	 * Query output dimension
	 */
	virtual int getOutputDim() = 0;
	
	/**
	 *
	 * Additional methods that are only for batch gradient algs
	 *
	 */
	
	/**
	 * Assuming we have accumulated a gradient somewhere, step in 
	 * that direction.
	 */
	virtual void batchStep() { NYI };
	
	/**
     * Calculate and accumulate instant gradient into a batch gradient
     * estimate
     */
	virtual void accumulateGrad(double, Observation&) { NYI };
	
	/**
	 * Reset the accumulated gradient only.
	 */
	virtual void resetGrad() { NYI };
	
	/**
	 * Update the search direction for the parameters from the current
	 * gradient and possible past information.
	 * For GPOMDP and most codes this is empty.
	 * This could hold conjugation code
	 * Or it can extract natural actor critic search direction
	 */
	virtual void computeDirection(int) { NYI };
	
	
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
