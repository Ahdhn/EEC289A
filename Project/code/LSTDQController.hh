#ifndef LSTDQController_hh
#define LSTDQController_hh

#include "PGBasics.hh"
#include "ValueController.hh"

namespace libpg {
    /**
     * This is a subclass of BasicController, implementing an abstract 
     * value controller. 
     */

    class LSTDQController : public ValueController {


    protected:

	// Value based method variables

	Policy* policy;      // Policy used to pick next action
	double tdDiscount;
	int actions;
    
	// LSTDQ specific variables

	int lstdqDimension;    // Size of basis vectors plus grads

	Vector previousPhi; // Copy of last basis feature (which serves as current phi)
	Vector currentPhi; // Copy of last basis feature (which serves as current phi)
	double lastReward;   // Copy of last reward
	double oldReward;    // The next to last reward
	bool instantStepRequired; // Are we waiting for the next obs to complete step?

	double discount; // Local copy of this. Always the same as controller discount

	Vector values; // Values computed by the approximator
	Matrix AInverse;
	Vector lstdqTrace;     // Trace variable for features
	Vector lstdqGrad;      // accumulated reward times obsTrace 
	Vector lstdqWeights; // computed search direction
	Vector lstdqTD;        // Temporal difference in basis, inc. grads
	Vector bVec;    // Search direction proposed by NAC.
	Matrix I;            // identity matrix shortcut

	int lstdqSteps; // The number of steps over which A and B have been computed

    public:
   
	LSTDQController(Approximator* approx, Policy* policy, double tdDiscount);
	virtual ~LSTDQController() {};

	virtual void discountTrace();

	virtual void getAction(Observation& obs, Vector& action, bool learn = true);
	virtual void getMostProbAction(Observation& obs, Vector& action) { NYI }
    
	virtual void resetTrace();
	virtual void resetGrad();
	virtual void accumulateGrad(Vector& reward, Observation& newObs);
	virtual void instantStep(Vector& reward);
	virtual void delayedInstantStep(double reward, int steps);
	virtual void computeDirection(int steps);
	virtual void batchStep();
	virtual void updateStatistics(double reward, int steps);
	virtual void scatter(Vector& v, Approximator::StatsEnum s);
	virtual void reduce(Vector& v, Approximator::StatsEnum s);
	virtual void setDiscount(double discount);
	void doApprox(Observation& obs, Vector& output) {approx->doApprox(obs, output);};

    };
}
#endif
