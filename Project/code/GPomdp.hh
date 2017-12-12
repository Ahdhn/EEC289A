#ifndef GPomdp_hh
#define GPomdp_hh
// $Id: GPomdp.hh 127 2007-09-10 17:20:04Z daa $

namespace libpg {
class GPomdp : public RLAlg {

public:

    /**
     * Don't call this one directly. Tears will result.
     */
    GPomdp() {}; 


    GPomdp(Controller* controller, Simulator* simulator, double discount, double stepSize);

    virtual double learnCore(int stepsPerEpoch, int& totalSteps);
    virtual double doSteps(Vector& totalReward, int steps, bool learn);

};
}
#endif
