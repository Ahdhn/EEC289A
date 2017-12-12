#ifndef RegressionController_hh
#define RegressionController_hh
// $Id: RegressionController.hh 60 2007-02-05 12:25:25Z daa $

namespace libpg {

/**
 * Instead of doing RL, it assumes the reward is the derivative of a loss function for the previous
 * observation (input) and action (approximator output) pair. Eligibility trace is cleared at every step.
 * Loss function derivative is computed in the Simulator::getReward() function.
 * Can be used in online or batch style
 */
class RegressionController : public Controller {

protected:

    Approximator* approx;
    Observation lastObs;
    Vector lastRewards;

public:

    RegressionController(Approximator* approx);
    virtual ~RegressionController() {};

    /**
     * Where the real coding is done to implement 
     */
    virtual void getAction(Observation& obs, Vector& action, bool computeGrad = true);
    virtual void computePartialsAndFeedback(Observation& obs, int sampledAct, double weight = 1.0);

    virtual void setLastReward(Vector& lastRewards);

    virtual void getMostProbAction(Observation& obs, Vector& action); 

    virtual int getInputDim();

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
    virtual int getNumParams();
    virtual void write(std::ostream& o);
    virtual void read(std::istream& o);
    virtual void scatter(Vector& v, Approximator::StatsEnum s);
    virtual void reduce(Vector& v, Approximator::StatsEnum s);

};
}
#endif
