#ifndef eGreedyPolicy_hh
#define eGreedyPolicy_hh

namespace libpg {
    /**
     * This is a subclass of Policy, implementing the e-Greedy
     * action selection.
     */

    
    class eGreedyPolicy : public Policy {

    public:

	/**
	 * The user should subclasse this abstract class in
	 * order to create his own epsilon function.
	 */
	class EpsilonFunction {
	    
	public:
	    
	    EpsilonFunction() {};
	    
	    /**
	     * Shut up compiler warning
	     */
	    virtual ~EpsilonFunction() {};
	    
	    /**
	     * Implements the epsilon function.
	     * @param steps number of steps done so far.
	     * @return temperature value.
	     */
	    virtual double getValue(int steps) = 0;
	    
	};
	
	
    private:
	
	double e;
	EpsilonFunction* epsilonFunc;
	
    public:


	eGreedyPolicy(double e = 0.0);
	eGreedyPolicy(EpsilonFunction* epsilonFunc);
	virtual ~eGreedyPolicy();

	/**
	 * Sample an action according to e-greedy, given distribution
	 * vector and observation.
	 */
	void getAction(Vector dist, Vector& action, Observation& obs);
    
    };
}
#endif
