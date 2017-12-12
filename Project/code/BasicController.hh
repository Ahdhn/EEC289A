#ifndef BasicController_hh
#define BasicController_hh

// $Id: BasicController.hh 127 2007-09-10 17:20:04Z daa $

namespace libpg {
/**
 * The first concrete controller class that only assumes that there is
 * exactly one ParameterisedDistribution that it talks to. Has not yet
 * made concrete the process of sampling from approximated
 * distributions.
 */

class BasicController : public Controller {


public:
   
    /**
     * the parameterised distribution
     */
    Approximator* approx;

    /**
     * The output of the approximator that we will pass through softmax
     */
    Vector dist;
    
    Vector distribution;

    BasicController(Approximator* approx);
    virtual ~BasicController() {};

    /**
     * Where the real coding is done to implement 
     */
    virtual void getAction(Observation& obs, Vector& action, bool computeGrad = true);
    virtual void getMostProbAction(Observation& obs, Vector& action);
    virtual void computePartialsAndFeedback(Observation& obs, int sampledAct, double weight = 1.0);

    virtual int  getInputDim();
    virtual int getNumParams();
    virtual void discountTrace();
    virtual void setDiscount(double discount);
    virtual void accumulateGrad(Vector& reward, Observation& newObs);
    virtual void resetGrad();
    virtual void instantStep(Vector& reward);
    virtual void batchStep();
    virtual void computeDirection(int steps);
    virtual void setStepSize(double stepSize);
    virtual void resetTrace();
    virtual void resetParams();
    virtual void randomizeParams(double maxRand);
    virtual double getMaxParam();
    virtual void write(std::ostream& o);
    virtual void read(std::istream& o);

    virtual void reduce(Vector& v, Approximator::StatsEnum s);
    virtual void scatter(Vector& v, Approximator::StatsEnum s);

    
};
}
#endif
