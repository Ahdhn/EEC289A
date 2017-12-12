#ifndef NeuralNetAtlas_hh
#define NeuralNetAtlas_hh
// $Id: NeuralNetAtlas.hh 127 2007-09-10 17:20:04Z daa $


namespace libpg {
#if (HAVE_BOOST_SANDBOX && HAVE_ATLAS)

class NeuralNetAtlas : public NeuralNetBatch {

protected:

public:

    NeuralNetAtlas(Vector& dims, Vector& squash) : NeuralNetBatch(dims, squash) {};
    virtual ~NeuralNetAtlas() {};

    virtual void doApprox(Observation& obs, Vector& output);
    virtual void feedbackGrad(Observation& obs, Vector& outputDeltas);
};

#endif
}
#endif
