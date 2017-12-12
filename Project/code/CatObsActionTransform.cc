/**
 * $Id: CatObsActionTransform.cc 127 2007-09-10 17:20:04Z daa $
 */

#include"PGBasics.hh"
#include"CatObsActionTransform.hh"

using namespace std;
namespace libpg {

CatObsActionTransform::CatObsActionTransform(Controller* controller,
					     size_t numRows, // = obs_dim,
					     size_t numCols,
					     size_t agents,
					     size_t actionDim,
					     size_t actions
					     ) : 
    TransformController(controller) {

    // Only one column is expected.
    assert(numCols==1);

    this->numRows = numRows;
    this->numCols = numCols;
    this->agents = agents;
    this->actionDim = actionDim;
    
    obsAction.init(numRows+actionDim,
		   numCols,
		   agents,
		   actions
		   );

    prevAction.resize(actionDim);
    prevAction.clear();
}


void CatObsActionTransform::getAction(Observation& obs, Vector& action, bool computeGrad) {
    
    // This observation window code can't handle multi-agent obs.
    assert(obs.getAgent() == 0);
    obsAction = obs;

    // Copy observations
    project( obsAction.getFeatures(),
	     range(0, numRows),
	     range(0, 1) )
	= project( obs.getFeatures(),
		   range(0, numRows),
		   range(0, 1) );

    // Append most recent obsAction.
    for (size_t i=0; i<actionDim; i++)
	(obsAction.getFeatures())(numRows+i, 0) = floor(prevAction(i));
    //cerr<<"CatObsAction="<<obsAction.features<<endl;

    try {
	controller->getAction(obsAction, action, computeGrad);
    } catch (std::exception& e) {
	cout<<"obsAction: "<<obsAction.getFeatures()<<endl;
	throw e;
    }

    prevAction.assign(action);

    //cerr<<"act="<<action<<endl;
    
}
}
