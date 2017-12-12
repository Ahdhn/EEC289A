#ifndef Bias_hh
#define Bias_hh
// $Id: BiasedApproximator.hh 87 2007-04-02 16:38:13Z buffet+daa $

#include "PGBasics.hh"
#include "TransformApproximator.hh"
#include "Teacher.hh"

namespace libpg {
    
    /**
     * This controller is a biased version of the ManyToOneController.
     */
    
    class Bias : public  TransformApproximator {
	
    protected: 
	
	Teacher* teacher;
	Vector teacherOutput;
	
    public:
	
	Bias(Approximator* approx, Teacher* teacher);
	virtual ~Bias() {};
	
	/**
	 * Where the real coding is done to implement 
	 */
	virtual void doApprox(Observation& obs, Vector& output);
		
    };
}
#endif
