#ifndef ImportanceSamplingController_hh
#define ImportanceSamplingController_hh

#include"BasicController.hh"

// $Id: ImportanceSamplingController.hh 29 2007-01-17 11:37:00Z buffet $

namespace libpg {

/**
 * A controller to guide the learning process with a "teacher" for a
 * period of time.
 *
 * It implements importance sampling to take into account the bias due
 * to the teacher.
 */
class ImportanceSamplingController : public TransformController {

public:

    enum ISType {noIS, weightedLocalIS, localOnTraceIS};

protected:

    ISType whichIS;

    double K; // Importance coefficient
    double Ka; // [olivier] let's say this is here for computing K's average
    int kSamp;
    double tRate; // Teaching rate (probability to use the teacher
                  // rather than the controller).
    bool isTeaching; // Are we in a teaching phase ?
    Teacher* teacher;
    BasicController* teacherController; // Turn teacher into a distribution

    double discount; // IS requires keeping a copy of the discount factor.

    Vector distribution;
    Vector teacherDistribution;
    bool needDistribution; // If using importance sampling, this tells
			   // that the probability distribution over
			   // actions should be stored in variable
			   // "distribution".



public:
   
    ImportanceSamplingController(Controller* controller, Teacher* teacher, int outputDim);
    virtual ~ImportanceSamplingController();

    virtual void getAction(Observation& obs, Vector& action, bool computeGrad = true);
    virtual void getMostProbAction(Observation& obs, Vector& action);
    virtual void computePartialsAndFeedback(Observation& obs, int sampledAct, double weight = 1.0) { NYI };
    
    virtual void instantStep(Vector& reward);
    /*virtual void batchStep();*/
    virtual void setDiscount(double discount);
    virtual void resetTrace();

    virtual Vector* getDistribution();

    void turnOffTeaching();
    void setTRate(double d);
    double getTRate();
    void setKa(double d);
    double getKa();
};
}
#endif
