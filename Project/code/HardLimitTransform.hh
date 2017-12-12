#ifndef HardLimitTransform_hh
#define HardLimitTransform_hh
// $Id: RegularizeTransform.hh 29 2006-12-19 02:16:20Z daa $

namespace libpg {

/**
 * Stop the underlying parameters (or values) from growing past a fixed absolute limit
 * Useful if you have an RL process that drives parameters past the point where exponentiating
 * them causes overflow problems.
 */
class HardLimitTransform : public TransformController {


protected:
   
    Controller* controller;
    double maxp;
    Vector params;
    
    virtual void checkParams();
public:

    HardLimitTransform(Controller* c, double penalty);
    ~HardLimitTransform() {};

    virtual void setMax(double maxp);
    virtual void instantStep(Vector& rewards);
    virtual void batchStep();

};
}
#endif
