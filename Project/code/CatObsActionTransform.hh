#ifndef CatObsActionTransform_hh
#define CatObsActionTransform_hh
/**
 * $Id
 */
namespace libpg {

/**
 * Concatenate current observation and previous action to produce a new
 * observation. Allows actions to be conditioned on previous action as well
 * as the observation.
 */
class CatObsActionTransform : public TransformController {


private:

    Observation obsAction; // Concatenated observation

    /**
     * Take care: values indicated below are for the observation
     * BEFORE transformation.
     */
    size_t numRows; // Number of observation dimensions
    size_t numCols; // Number of observation columns
    size_t agents; // Number of agents
    size_t actionDim; // Number of action dimensions
    
    Vector prevAction; // The action that generated the current obs

public:
    
    /** 
     * All of the parameter information here is availalble from the simulator
     * class.
     */
    CatObsActionTransform(Controller* controller, 
			  size_t numRows, // = obs_dim,
			  size_t numCols, 
			  size_t agents,
			  size_t actionDim,
			  size_t actions);

    ~CatObsActionTransform() {};

    virtual void getAction(Observation& obs, 
			   Vector& action, 
			   bool computeGrad = true);
};
}
#endif
