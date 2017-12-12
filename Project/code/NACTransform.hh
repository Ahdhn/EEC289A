#ifndef NACTransform_hh
#define NACTransform_hh

// $Id: NACTransform.hh 127 2007-09-10 17:20:04Z daa $

namespace libpg {

    class NACTransform : public TransformController {


    protected:

	int nacDimension;    // Size of basis vectors plus grads

	Observation lastObs; // Copy of last observation
	double lastReward;   // Copy of last reward
	bool instantStepRequired; // Are we waiting for the next obs to complete step?

	double discount; // Local copy of this. Always the same as controller discount

	Matrix AInverse;
	Vector nacTrace;     // Trace variable for features
	Vector nacGrad;      // accumulated reward times obsTrace 
	Vector nacDirection; // computed search direction
	Vector nacTD;        // Temporal difference in basis, inc. grads
	Vector nacParamUpdate;    // Search direction proposed by NAC.
	Matrix I;            // identity matrix shortcut
    
	double tdDiscount;   // Critic discount factor
    
    public:
    
	NACTransform(Controller* controller, double tdDiscount);


	virtual void discountTrace();

	virtual void getAction(Observation& obs, Vector& action, bool computeGrad = true);
	virtual void getMostProbAction(Observation& obs, Vector& action) { NYI }

	virtual void resetTrace();
	virtual void resetGrad();
	virtual void accumulateGrad(double reward, Observation& newObs);
	virtual void instantStep(Vector& reward);
	virtual void delayedInstantStep(double reward, Observation& newObs);
	virtual void computeDirection(int steps);
	virtual void batchStep();
	virtual void updateStatistics(double reward, Observation& newObs);
	virtual void scatter(Vector& v, Approximator::StatsEnum s);
	virtual void reduce(Vector& v, Approximator::StatsEnum s);
	virtual void setDiscount(double discount);
    

    };
}
#endif
