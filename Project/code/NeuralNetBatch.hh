#ifndef NeuralNetBatch_hh
#define NeuralNetBatch_hh
// $Id: NeuralNetBatch.hh 127 2007-09-10 17:20:04Z daa $


#include"NeuralNet.hh"

namespace libpg {
/**
 * A version of neural network that supports aggregated gradients
 */

class NeuralNetBatch : public NeuralNet {

public:

    Matrix* grads;

    
    NeuralNetBatch(Vector& dims, Vector& squash);
    NeuralNetBatch(int in, int out);

    virtual void batchStep();
    virtual void computeDirection(int steps);
    virtual void accumulateGrad(double reward, Observation& newObs);
    virtual void resetGrad();


    virtual void reduce(Vector& v, StatsEnum s);
    virtual void scatter(Vector& v, StatsEnum s);

};
}
#endif
