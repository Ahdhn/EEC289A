/**
 * $Id: ObsWindowTransform.cc 127 2007-09-10 17:20:04Z daa $
 */

#include"PGBasics.hh"
#include"ObsWindowTransform.hh"

namespace libpg {

ObsWindowTransform::ObsWindowTransform(Controller* controller,
			   size_t maxHistory,
			   size_t numRows, // = obs_dim,
			   size_t numCols,
			   size_t agents
			   ) : 
    TransformController(controller) {


    this->maxHistory = maxHistory;
    this->numRows = numRows;
    this->numCols = numCols;
    this->agents = agents;
    
    history.init(numRows*maxHistory,
		 numCols,
		 agents,
		 1  // History doesn't need to know the action dimension
		 );
    history.getFeatures().clear();
}


void ObsWindowTransform::getAction(Observation& obs, Vector& action, bool computeGrad) {
    
    // This observation window code can't handle multi-agent obs.
    assert(obs.getAgent() == 0);
    history.setAgent(obs.getAgent());
    history.getEligible() = obs.getEligible();

    // Translate the history
    project( history.getFeatures(),
	     range(0, numRows*(maxHistory-1)),
	     range(0,numCols) )
	= project( history.getFeatures(),
		   range(numRows, numRows*maxHistory),
		   range(0 ,numCols) );

    // Append most recent history.
    project( history.getFeatures(),
	     range(numRows*(maxHistory-1), numRows*maxHistory),
	     range(0,numCols) )
	= obs.getFeatures();
    
    /*cerr<<" "<<endl;
      cerr<<obs.features<<endl;*/
    //cerr<<"obsWindow="<<history.features<<endl;
    
    try {
	controller->getAction(history, action, computeGrad);
    } catch (std::exception& e) {
	//cout<<"history: "<<history.features<<endl;
	throw e;
    }

}
}
