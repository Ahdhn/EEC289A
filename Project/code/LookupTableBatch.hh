#ifndef LookupTableBatch_hh
#define LookupTableBatch_hh
// $Id: LookupTableBatch.hh 127 2007-09-10 17:20:04Z daa $

#include"LookupTable.hh"

namespace libpg {
/**
 * A version of neural network that supports aggregated gradients
 */

class LookupTableBatch : public LookupTable {

protected:

    Matrix grads;
    LookupTableBatch() {};

public:

    LookupTableBatch(int observations, int outputs);
    virtual ~LookupTableBatch();
    
    virtual void batchStep();
    virtual void computeDirection(int steps);
    virtual void accumulateGrad(double reward, Observation& newObs);
    virtual void resetGrad();

    virtual void scatter(Vector& v, Approximator::StatsEnum s);
    virtual void reduce(Vector& v, Approximator::StatsEnum s);
};
}
#endif
