#ifndef Teacher_hh
#define Teacher_hh
// $Id: Teacher.hh 0 2007-01-17 10:59:00Z buffet $

namespace libpg {


/**
 * Decision-maker used as a teacher for an RL algorithm.
 */



/**
 * Required by importance sampling policy-gradients and Biased Controllers
 */

    class Teacher : public Approximator {
	
    public:
	
	/**
	 * Shut up compiler warning
	 */
	virtual ~Teacher() {};
		
	/**
	 * Return a set of likelihoods over actions. May be negative.
	 * Imagine these values are being added directly to the learning approximators values.
	 */
	virtual void doApprox(Observation& obs, Vector& output) = 0;

	/**
	 * Need to implement these two on a case by case basis
	 */
	virtual int getInputDim() = 0;
	virtual int getOutputDim() = 0;

	
	/**
	 * Teacher has no parameters to learn by definition.
	 */
	virtual int getNumParams() { return 0; };

	/**
	 * No learning, so can't feedback.Generate error if called.
	 */
	virtual void feedbackGrad(Observation& obs, Vector& deltas) { NYI };

	/**
	 * Again, no learning going on, so this should generate an error.
	 */
	virtual void discountTrace() { NYI };

	/**
	 * Again, no learning going on, so this should generate an error.
	 */
	virtual void setDiscount(double) { NYI };


	/**
	 * Again, no learning going on, so this should generate an error.
	 */
	virtual void instantStep(double) { NYI };


	/**
	 * Again, no learning going on, so this should generate an error.
	 */
	virtual void setStepSize(double) { NYI };

	/**
	 * Again, no learning going on, so this should generate an error.
	 */
	virtual void resetTrace() { NYI };
	virtual void resetParams() { NYI };


	/**
	 * Again, no learning going on, so this should generate an error.
	 */
	virtual void randomizeParams(double) { NYI };

	virtual double getMaxParam() { return 0.0; }


    };
}
#endif
