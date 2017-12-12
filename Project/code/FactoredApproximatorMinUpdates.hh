#ifndef FactoredApproximatorMinUpdates_hh
#define FactoredApproximatorMinUpdates_hh
// $Id: FactoredApproximatorMinUpdates.hh 29 2006-12-19 02:16:20Z daa $

#include"FactoredApproximator.hh"

namespace libpg {

/**
 * An approximator that prevents uncessesary 0 updates in a many to one situation.
 */

class FactoredApproximatorMinUpdates : public FactoredApproximator {


public:

    double discount;

    bVector* eligible;
    int currentTime;
    Vector lastUpdateTime;
    bVector hasBeenEligible;
    bVector hasBeenEligibleLongTerm;

    FactoredApproximatorMinUpdates(Approximators approximators, bool splitObs);

    /**
     * Where the real coding is done to implement 
     */
    virtual void doApprox(Observation& obs, Vector& outputs);

    virtual void discountTrace();
    virtual void setDiscount(double discount);
    virtual void instantStep(double reward);
    virtual void resetTrace();

};
}
#endif
