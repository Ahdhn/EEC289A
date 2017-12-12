#ifndef TransformApproximator_hh
#define TransformApproximator_hh
// $Id: Approximator.hh 48 2007-02-01 06:52:57Z buffet $

#include "PGBasics.hh"

namespace libpg {

/**
 * Generic transformation of an approximator to permit screwing with the output values
 * of approximators without having to reimplement every single method of an approximator.
 */

/**
 * Knows about, and performs, basic policy, or value, gradient operations.
 */
    
    class TransformApproximator : public Approximator {
	
	
    protected: 
	Approximator* approx;
	
    public:
	
	
	/**
	 * Create a wrapper around an existing approximator only.
	 * @param  approx base approximator
	 */
	TransformApproximator(Approximator* approx);

	/**
	 * Shut up compiler warning
	 */
	virtual ~TransformApproximator() {};
	
	/** 
	 * Every single method here just calls its base equivalent.
	 */
	virtual int getNumParams();
	virtual void doApprox(Observation& obs, Vector& output);
	virtual void feedbackGrad(Observation& obs, Vector& deltas); 
	virtual void discountTrace();
	virtual void setDiscount(double);
	virtual void instantStep(double);
	virtual void setStepSize(double);
	virtual void resetTrace();
	virtual void resetParams();
	virtual void randomizeParams(double);
	virtual double getMaxParam();
	virtual void write(std::ostream&);
	virtual void read(std::istream&);
	virtual int getInputDim();
	virtual int getOutputDim();
	virtual void batchStep();
	virtual void accumulateGrad(double, Observation&);
	virtual void resetGrad();
	virtual void computeDirection(int);
	virtual void reduce(Vector&, Approximator::StatsEnum);
	virtual void scatter(Vector&, Approximator::StatsEnum);
	
    };
}
#endif
