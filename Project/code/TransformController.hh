#ifndef TransformController_hh
#define TransformController_hh

// $Id: TransformController.hh 127 2007-09-10 17:20:04Z daa $

namespace libpg {

/**
 * A dummy controller that serves as a basis for sub-classes to
 * build transformations of observations and actions on without having
 * to define all the abstract methods below.
 */
class TransformController : public Controller {


protected:

    Controller* controller;


public:
   
    TransformController(Controller* controller);
    virtual ~TransformController() {};

    virtual void getAction(Observation& obs, Vector& action, bool computeGrad = true);
    virtual void computePartialsAndFeedback(Observation& obs, int sampledAct, double weight = 1.0);
    virtual void getMostProbAction(Observation& obs, Vector& action);
    
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
