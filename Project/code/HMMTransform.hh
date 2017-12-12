#ifndef HMMTransform_hh
#define HMMTransform_hh
/**
 * $Id
 */

namespace libpg {
/**
 * Take observations and turn them into HMM prediction vector.
 * Optionally learn HMM online.
 */

#define TERM_PROB_THRESH 0.01
#define MIN_PROB         0.001

class HMMTransform : public TransformController {


private:

    boost::numeric::ublas::vector<int> oHistory; // observation history
    boost::numeric::ublas::vector<int> aHistory; // action history
    size_t maxHistory; // Length of history to work with.
    size_t window; // Length of the window used to compute belief state
    size_t historyLength; // How long is the history right now ?
    size_t numStates; // Assumed number of states
    size_t numActions; // Number of action symbols
    size_t numObs; // Number of observation symbols
    bool learn; // Are we learning to track the HMM
    
    size_t prevAction; // The action that generated the current obs

    Observation hmm;

    Matrix alpha;
    Matrix beta;
    Matrix gamma;
    Matrix* xi;
    Matrix* A;
    Matrix* B;
    Vector pi;

    /**
     * Do parameter restimation
     */
    void reEstimate(void);
    void estimate(void);

    /**
     * Prob est routines
     */
    double forwardEst(Vector &scale, size_t L);
    void   backwardEst(Vector &scale, size_t L);

    /**
     * Counting transitions
     */
    void gammaEst(size_t L);
    void xiEst(void);

    /**
     * Maximisation
     */
    void maximiseTrans(void);
    void maximiseDiscreteObs(void);

public:
    
    HMMTransform(Controller* controller, 
		 size_t maxHistory, size_t window,
		 size_t numStates,
		 size_t numObs,
		 size_t numActions);

    ~HMMTransform();

    virtual void getAction(Observation& obs, 
			   Vector& action, 
			   bool computeGrad = true);
};
}
#endif
