#ifndef LookupTable_hh
#define LookupTable_hh
// $Id: LookupTable.hh 127 2007-09-10 17:20:04Z daa $

namespace libpg {

/**
 * Simplest possible parameterised distribution. Observation indexes a
 * column. Column returns real value that is likelihood of each action.
 */

class LookupTable : public Approximator {

protected:

    /**
     *Constructor for specialisations. 
     */
    LookupTable () { }

    /**
     * number of actions (output dimension)
     */
    int actions;
    
    /**
     * Number of observations
     */
    int observations;
    
    /**
     * parameter matrix
     */
    Matrix params;
    
    /**
     * Eligibility trace;
     */
    Matrix trace;
    
    /**
     * Discount factor
     */
    double discount;
    
    /**
     * Step size
     */
    double stepSize;
    
public:
   
    LookupTable(int observations, int outputs);
    virtual ~LookupTable() {};

    virtual void doApprox(Observation& obs, Vector& outputs);
    virtual void feedbackGrad(Observation& obs, Vector& deltas); 
    virtual void discountTrace();
    virtual void setDiscount(double discount);
    virtual void instantStep(double reward);
    virtual void setStepSize(double stepSize);
    virtual void resetTrace();
    virtual void resetParams();
    virtual void randomizeParams(double maxRand);
    virtual double getMaxParam();
    virtual void write(std::ostream& o);
    virtual void read(std::istream& o);
    virtual int getNumParams ();
    virtual int getInputDim();
    virtual int getOutputDim();

    virtual void scatter(Vector& v, Approximator::StatsEnum s);
    virtual void reduce(Vector& v, Approximator::StatsEnum s);
};
}
#endif
