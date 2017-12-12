#ifndef CyclicPolicyBias_hh
#define CyclicPolicyBias_hh
// $Id: CyclicPolicyBias.hh 127 2007-09-10 17:20:04Z daa $

namespace libpg {

/**
 * Knows about time using the clock from RLAlg in the features.
 * Adds a constant to each approximator output based on the clock,
 * and a phase time. Requires number of approximator outputs be 
 * a factor of the cycle time.
 */
class CyclicPolicyBias : public Approximator {

private:

    int parameters;

    Approximator* approx;
    int stepsPerControl; // How many consecutive steps to add bias to a particular control
    int cycleTime; // Total length of policy cycle
    double bias; // How much to add to outputs
    
public:
   
    /**
     * Constructor, do whatever you need 
     */
    CyclicPolicyBias(Approximator*, int cycleTime, double bias);

    /**
     * Adds constant to normal approximator
     */
    virtual void doApprox(Observation& obs, Vector& output);


    /**
     * Everything else is just a pass through to the normal approximator
     */
    virtual void feedbackGrad(Observation& obs, Vector& deltas); 
    virtual void discountTrace();
    virtual void setDiscount(double discount) ;
    virtual void instantStep(double reward);
    virtual void setStepSize(double stepSize);
    virtual void resetTrace();
    virtual void resetParams();
    virtual void randomizeParams(double maxRand);
    virtual double getMaxParam();
    virtual void write(std::ostream& o);
    virtual void read(std::istream& o);
    virtual void batchStep();
    virtual void accumulateGrad(double reward, Observation& newObs);
    virtual void resetGrad();
    virtual void computeDirection(int steps);
    virtual int getNumParams();

    virtual int getInputDim();
    virtual int getOutputDim();

    virtual void reduce(Vector& v, Approximator::StatsEnum s);
    virtual void scatter(Vector& v, Approximator::StatsEnum s);

};
}
#endif
