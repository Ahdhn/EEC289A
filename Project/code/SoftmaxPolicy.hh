#ifndef SoftmaxPolicy_hh
#define SoftmaxPolicy_hh

namespace libpg {
    /**
     * This is a subclass of Policy, implementing the Softmax
     * action selection, with Gibbs (or Boltzmann) distribution.
     */

    class SoftmaxPolicy : public Policy {

    public:


	/**
	 * The user should subclasse this abstract temperature
	 * class in order to create his own temperature decaying
	 * function.
	 */
	class Temperature {

	public:

	    Temperature() {};

	    /**
	     * Shut up compiler warning
	     */
	    virtual ~Temperature() {};

	    /**
	     * Implements the temperature decaying function.
	     * @param steps number of steps done so far.
	     * @return temperature value.
	     */
	    virtual double getValue(int steps) = 0;

	};


	SoftmaxPolicy(double constTemp = 1.0);
	SoftmaxPolicy(Temperature* tempFunc);
	virtual ~SoftmaxPolicy();

	/**
	 * Sample an action according to Softmax, given distribution
	 * vector and observation.
	 */
	void getAction(Vector dist, Vector& action, Observation& obs);


    private:

	double temp;
	Temperature* tempFunc;
    
    };
}
#endif
