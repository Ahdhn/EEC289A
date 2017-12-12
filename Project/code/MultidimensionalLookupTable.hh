#ifndef MultidimensionalLookupTable_hh
#define MultidimensionalLookupTable_hh

namespace libpg {

/**
 * A lookup table that accepts as inputs a finite windows of past history
 */
class MultidimensionalLookupTable : public LookupTableBatch {
    
protected:
    int maxHistory;
    int* powerArray;
    
public:

    MultidimensionalLookupTable (int observations, int outputs, int memoryLength);
    
    virtual void doApprox (Observation& obs, Vector& output);    
    virtual void feedbackGrad (Observation& obs, Vector& deltas);
    virtual int getInputDim();
    
};
}
#endif
