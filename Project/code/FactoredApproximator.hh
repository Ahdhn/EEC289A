// $Id$

/**
 * Factored Approximator to aggregate independent approximators into one.
 * Replaces old ManyToOneController concept.
 */

#ifndef FactoredApproximator_hh
#define FactoredApproximator_hh

#include"PGBasics.hh"

namespace libpg {


    class FactoredApproximator : public Approximator {
	
    public:
	typedef std::vector<Approximator*> Approximators;
	
    protected:
	Approximators approximators;
	Vector dummyDist;
	bool splitObs;

    public:

	
	FactoredApproximator(Approximators approximators, bool splitObs);
	
	/**
	 * Shut up compiler warning
	 */
	virtual ~FactoredApproximator() {};
       
	
	/**
	 * Constructor, do whatever you need to.
	 * Approximator() {};
	 */
	virtual int getNumParams();
	
	/**
	 * Aggregate the output of each approximator
	 */
	virtual void doApprox(Observation& obs, Vector& output);
	
	/**
	 * Feed back a gradient signal to each approximator. Final gradient accumulated directly to
	 * traces in each approximator.
	 */
	virtual void feedbackGrad(Observation& obs, Vector& deltas) ; 
	
	/**
	 * Just multiply trace by  discount factors
	 */
	virtual void discountTrace();
	
	/**
	 * Set the discount factors. Does not necessarily have to be constant.
	 */
	virtual void setDiscount(double);
	
	/**
	 * Calculate instant gradient and step the parameters
	 */
	virtual void instantStep(double);
	
	/**
	 * Set the step size. Does not necessarily have to be a constant.
	 */
	virtual void setStepSize(double);
	
	/**
	 * Reset the eligibility trace to 0
	 */
	virtual void resetTrace();
	
	/**
	 * Do whatever resetting the params means.
	 */
	virtual void resetParams();
	
	/**
	 * Create a new random set of params.
	 */
	virtual void randomizeParams(double);
	
	/**
	 * Retrieve largest absolute parameter value.
	 */
	virtual double getMaxParam();
	
	/**
	 * Read or write parameters
	 */
	virtual void write(std::ostream&);
	virtual void read(std::istream&);
	
	
	/**
	 * Query expected number of inputs
	 */
	virtual int getInputDim();
	
	/**
	 * Query output dimension
	 */
	virtual int getOutputDim();
	
	/**
	 * Additional methods that are only for batch gradient algs
	 */
	
	/**
	 * Assuming we have accumulated a gradient somewhere, step in 
	 * that direction.
	 */
	virtual void batchStep();
	
	/**
	 * Calculate and accumulate instant gradient into a batch gradient
	 * estimate
	 */
	virtual void accumulateGrad(double, Observation&);
	
	/**
	 * Reset the accumulated gradient only.
	 */
	virtual void resetGrad();
	
	/**
	 * Update the search direction for the parameters from the current
	 * gradient and possible past information.
	 * For GPOMDP and most codes this is empty.
	 * This could hold conjugation code
	 * Or it can extract natural actor critic search direction
	 */
	virtual void computeDirection(int);
	
	
	/**
	 * Add either the parameters or gradient or traces into a vector 
	 * @param v vector to store stuff into. v must be at least as big
	 * as total params/grads/traces. Will not be resized. Will not be cleared.
	 * @param s whether to collect grads or params or traces
	 */
	virtual void reduce(Vector&, Approximator::StatsEnum);
	
	
	/**
	 * Disperse either the parameters or gradient or any other
	 * type of stat into appropriate data structures for this
	 * class. Actually defers to sub approximators.
	 * @param v Vector to replace current params or grads
	 * or traces. Must be correct size 
	 * @param s whether to restore
	 * grads or params or traces
	 */
	virtual void scatter(Vector&, Approximator::StatsEnum);
	
    };


}

#endif
