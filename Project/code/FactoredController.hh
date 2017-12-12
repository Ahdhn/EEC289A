#ifndef FactoredController_hh
#define FactoredController_hh
// $Id: FactoredController.hh 127 2007-09-10 17:20:04Z daa $


namespace libpg {
/**
 * A controller that is really a collection of other controllers. This
 * controller simply fields observations to all controllers, and
 * gathers back a single dimension action from each of them.
 */

class FactoredController : public Controller {


public:

    /**
     * the sub controllers
     */
    typedef std::vector<Controller*> Controllers;
    Controllers controllers;

    /**
     * Does each controller get its own observation
     */
    bool splitObs;
    Vector dummyAction;
    Vector subReward;

    /**
     * Is there a reward for each controller?
     */
    bool localRewards;

    FactoredController(Controllers controllers, bool splitObs, bool localRewards);
    virtual ~FactoredController() {};

    /**
     * Where the real coding is done to implement 
     */
    virtual void getAction(Observation& obs, Vector& action, bool computeGrad = true);
    virtual void getMostProbAction(Observation& obs, Vector& action);

    virtual void computePartialsAndFeedback(Observation& obs, int sampledAct, double weight = 1.0) { NYI };
    virtual int getInputDim();
    virtual int getNumParams();
    virtual void discountTrace();
    virtual void setDiscount(double discount);
    virtual void accumulateGrad(Vector& rewards, Observation& newObs);
    virtual void resetGrad();
    virtual void instantStep(Vector& rewards);
    virtual void computeDirection(int steps);
    virtual void batchStep();
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
