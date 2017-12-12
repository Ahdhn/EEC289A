#ifndef OLPomdp_hh
#define OLPomdp_hh
// $Id: OLPomdp.hh 127 2007-09-10 17:20:04Z daa $

namespace libpg {
class OLPomdp : public RLAlg {

public:

    
    /**
     * Don't call this one. Tears will result.
     */
    OLPomdp() {}; 

    OLPomdp(Controller* controller, Simulator* simulator, double discount, double stepSize);
    virtual double doSteps(Vector& totalReward, int steps, bool learn);

};
}
#endif
