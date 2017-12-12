#ifndef RegularizeTransform_hh
#define RegularizeTransform_hh
// $Id: RegularizeTransform.hh 127 2007-09-10 17:20:04Z daa $

namespace libpg {
class RegularizeTransform : public TransformController {


protected:
    
    double stepSize;
    double penalty;
    Vector params;
    Vector penaltyVector;
    
    void penalise();

public:

    RegularizeTransform(Controller* c, double penalty);
    ~RegularizeTransform() {};

    virtual void setPenalty(double penalty);
    virtual void setStepSize(double stepSize);
    virtual void instantStep(Vector& rewards);
    virtual void batchStep();

};
}
#endif
