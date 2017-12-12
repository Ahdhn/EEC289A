#ifndef LineSearchAlg_hh
#define LineSearchAlg_hh
// $Id: LineSearchAlg.hh 127 2007-09-10 17:20:04Z daa $

namespace libpg {

/**
 * Perform a line search given a batch estimate coputed by GPomdp. 
 * Line search is controlled from the learnCore level, so the dosteps and learn functions are carried over nicely from GPomdp.
 */

class LineSearchAlg : public GPomdp {

public:

    int lineSearchSteps; // number of steps to evaluate potential step sizes for.
    int maxTries; // Maximum number of loops line search will go for
    double tolerance; // Quit when gradient magnitude falls below this
    double exponentBase; // Multiplicative factor for each step in exponential line search
    double backOff; // How far to back off at the start of each line search
    double stepSizeBackup; // Original stepsize
    double downhilllThresh; // How much to tolerate downhilll steps.

    Vector origParams; // Parameter vector at start of line search
    Vector params; // Curent parameter vector

    Vector grads; // Gradient returned from formal computation.
    Vector testGrads; // Approximate gradient returned after line search eval

    int seed;

    /**
     * Don't call this one directly. Tears will result.
     */
    LineSearchAlg() {}; 
    
    LineSearchAlg(Controller* controller, 
		  Simulator* simulator, 
		  double discount, 
		  double stepsize, 
		  int lineSearchSteps, 
		  double tolerance, 
		  double backOff, 
		  double exponentBase,
		  double downhilllThresh,
		  int maxTries);


    virtual void lineSearchInit(int lineSearchSteps, 
				double tolerance, 
				double backOff, 
				double exponentBase,
				double downhilllThresh,
				int maxTries);
    
    virtual double learnCore(int stepsPerEpoch, int& totalSteps);
    virtual double lineSearch(int& totalSteps, double lastEpochVal, Vector& totalRewards);
    virtual double incrementLineSearch(double currentStep, Vector& totalRewards, int& totalSteps, double& currentReward);


};
}
#endif
