#ifndef NeuralNet_hh
#define NeuralNet_hh
// $Id: NeuralNet.hh 127 2007-09-10 17:20:04Z daa $

namespace libpg {
class NeuralNet : public Approximator {

protected:

    int layers; // Number of layers of parameters
    int parameters; // Number of parameters

    Vector squash; // Which layer outputs to squash
    Vector dims; // Size of each layer.

    Matrix* layerParams; // The parameters of the neural network
    Matrix* layerTraces; // The trace vectors
    Vector* activations; // activations (outputs) of each layer except output and input
    Vector* deltas;      // Fed back errors at each layer except output and input

    double discount; // Discount factor
    double stepSize; // step size for online steps

public:

    NeuralNet(int inputs, int outputs);
    NeuralNet(Vector& dims, Vector& squash);
    virtual ~NeuralNet();

    virtual void init();

    virtual void doApprox(Observation& obs, Vector& output);
    virtual void feedbackGrad(Observation& obs, Vector& outputDeltas); 
    
    void squashVec(Vector& x);
    void dSquashVec(Vector& x, Vector& v);

    virtual void discountTrace();
    virtual void setDiscount(double discount);
    virtual void instantStep(double reward);
    virtual void setStepSize(double stepSize);
    virtual void resetTrace();
    virtual void resetParams();
    virtual void randomizeParams(double maxRand);
    virtual double getMaxParam();
    virtual int getNumParams();
    virtual void write(std::ostream& o);
    virtual void read(std::istream& o);
    virtual int getInputDim();
    virtual int getOutputDim();

    virtual void reduce(Vector& v, Approximator::StatsEnum s);
    virtual void scatter(Vector& v, Approximator::StatsEnum s);
};
}
#endif
