#ifndef QLearningController_hh
#define QLearningController_hh

namespace libpg {
/**
 * This is a subclass of ValueController, implementing a Watkin's
 * Q-Learning Labda specific controller. 
 */

class QLearningController : public ValueController {


private:

    Policy* policy;
    int lastAct;         // Last action taken
    Observation lastObs; // Last observation
    double error;        // Temporal difference error
    double lambda;     // Discount (lambda)
    double tdDiscount;   // Temporal difference discount
    double lastRew;      // Previous reward
    bool computeGrad;    // Should we learn by experience?


    Policy* greedyPolicy; // Greedy Policy
    Vector greedyAct;     // Greedy action
    Vector dist;
    Vector lastDist;      // Will hold previous values
    

public:
   
    QLearningController(Approximator* approx, Policy* policy, double tdDiscount);
    virtual ~QLearningController();

    /**
     * Where the real coding is done to implement Q learning error calculation
     */
    virtual void getAction(Observation& obs, Vector& action, bool computeGrad = true);

    virtual void instantStep(Vector& reward);

    virtual void discountTrace();


    virtual void setDiscount(double lambda);


};
}
#endif
