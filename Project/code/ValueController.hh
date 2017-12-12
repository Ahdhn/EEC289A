#ifndef ValueController_hh
#define ValueController_hh

#include "PGBasics.hh"
#include "BasicController.hh"

namespace libpg {
    /**
     * This is a subclass of BasicController, implementing an abstract 
     * value controller. 
     */

    class ValueController : public BasicController {


    public:
	
	ValueController(Approximator* approx) : BasicController(approx) {};

	virtual void getAction(Observation& obs, Vector& action, bool computeGrad = true) = 0;
	virtual void getMostProbAction(Observation& obs, Vector& action) { NYI }

	// No gradients to compute.
	virtual void computePartialsAndFeedback(Observation& obs, int sampledAct, double weight = 1.0) { NYI };

    };
}
#endif
