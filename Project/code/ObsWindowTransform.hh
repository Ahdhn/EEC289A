#ifndef ObsWindowTransform_hh
#define ObsWindowTransform_hh
/**
 * $Id
 */

namespace libpg {
/**
 * Take observations and stack them in an history vector of fixed length.
 * History vector then becomes the observation vector.
 */



class ObsWindowTransform : public TransformController {


private:

    Observation history; // Current history
    size_t maxHistory; // Length of history to work with.

    /**
     * Take care: values indicated below are for the observation
     * BEFORE transformation.
     */
    size_t numRows; // Number of observation dimensions
    size_t numCols; // Number of observation columns
    size_t agents; // Number of agents
    
    Vector prevAction; // The action that generated the current obs

    void updateHistory(size_t row, size_t action, size_t obsx);

public:
    
    ObsWindowTransform(Controller* controller, 
		       size_t maxHistory,
		       size_t numRows, // = obs_dim,
		       size_t numCols, 
		       size_t agents);

    ~ObsWindowTransform() {};

    virtual void getAction(Observation& obs, 
			   Vector& action, 
			   bool computeGrad = true);
};
}
#endif
