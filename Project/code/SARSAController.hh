#ifndef SARSAController_hh
#define SARSAController_hh

namespace libpg {
/**
 * This is a subclass of ValueController, implementing a SARSA Lambda
 * specific controller. 
 */

class SARSAController : public ValueController {

protected:

    Policy* policy;
    int lastAct;         // Last action taken
    Observation lastObs; // Last observation
    double error;        // Temporal difference error
    double lambda;     // Discount (lambda)
    double tdDiscount;   // Temporal difference discount
    double lastRew;      // Previous reward
    bool computeGrad;    // Should we learn by experience?

    Vector dist;
    Vector lastDist;     // Holds old and new transition values

public:
   
    SARSAController(Approximator* approx, Policy* policy, double tdDiscount);
    virtual ~SARSAController();


    /**
     * Where the real coding is done to implement the TD error calculation
     */
    virtual void getAction(Observation& obs, Vector& action, bool computeGrad = true);

    virtual void instantStep(Vector& reward);

    virtual void discountTrace();
    virtual void setDiscount(double lambda);

};
}
#endif
