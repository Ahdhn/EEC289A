/**
 * $Id: PSRTransform.cc 127 2007-09-10 17:20:04Z daa $
 */

#include"PGBasics.hh"
#include"PSRTransform.hh"

using namespace std;

namespace libpg {
/**
 * These defines are not tunable parameters.
 */
#define ACCURACY_DISCOUNT 0.99
#define DIGITS 10
#define ADD_ALL_CORES false
#define PSR_CONF_FREQ 10000
#define PSR_BINARY_THRESH 0.5
#define APPEND_OBS false
#define BINARY_OBS false

#define ADDCHRS(str, a, o) {str.push_back((char)a + '0'); \
	str.push_back((char)o + '0');}

#define PREVROW(row) ((row + maxHistory - 1)%maxHistory)

PSRTransform::PSRTransform(Controller* controller, 
			   size_t maxHistory, 
			   size_t numObs, 
			   size_t numActions,
			   double stepSize,
			   int passesBeforeAdd,
			   double minKappaForCore) : 

    TransformController(controller) {

    this->maxHistory = maxHistory;
    this->numObs = numObs;
    this->numActions = numActions;
    this->stepSize = stepSize;
    this->passesBeforeAdd = passesBeforeAdd;
    this->minKappaForCore = minKappaForCore;

    appendObs = APPEND_OBS;

    // because strings are currently used to represent history we
    // can't have more than 10 observations or actions, since
    // assumption is made that one character is one observation or
    // action. In the future we could have a specialised vector to
    // implement the buffer. Strings are easy because we have instant
    // substring and comparison operations.
    assert(numObs <= DIGITS);
    assert(numActions <= DIGITS);

    tests[nullTest] = 0;
    // This will also create appropriate size testPrediction matrix
    maxTestLen=0;
    testPredictions.resize(maxHistory, 1);
    testPredictions.clear();
    testPredictions(0,0) = 1.0;
    TestList tl;
    tl.push_back(tests.find(nullTest));
    completeTestSet(0, tl);

    // coreTest set starts with just null test
    coreTests[nullTest] = 0;

    params.resize(coreTests.size());
    params.clear();
    Q.resize(maxHistory, coreTests.size());
    Q.clear();
    testHistory.resize(maxHistory);
    testHistory.clear();
    psrVector.resize(coreTests.size()); psrVector.clear();

    // PSR vector size must be fixed to input dim of controller.
    // Note that PSR type is Observation.
    psrObs.init(controller->getInputDim(), 1, 1, 1);
    maxCores = controller->getInputDim() - (appendObs?numObs:0);

  
    discover = true;
    learn = true;
    addAllCores = ADD_ALL_CORES;
    binaryObs = BINARY_OBS;
    predictionAccuracy = 0.0;
    confidence = 0.0;
    lastRow=0;
    prevAction = 0;
    passes=0;

}


/**
 * Add extended tests given a set of new core tests. 
 * The new core tests are already memebers of the test set
 * but we need to add all a,o prefixes of all core tests,
 * and we also need to ensure that the length-2 substring test
 * is also a test, and that all obs/action suffix extensions
 * are also a test.
 */
void PSRTransform::completeTestSet(size_t row, TestList& newCores) {

    TestsMap::iterator test;
    std::string candidate;
    size_t initSize = tests.size();
    int pos;
    bool added;

    // Pass over all new cores and add the prefix
    // versions if they are not already there.
    for (TestList::iterator nc=newCores.begin(); nc != newCores.end(); nc++) {

	// Check for existence of prefix versions. Add if not found.
	for (char a=0; a < (char)numActions; a++) {
	    for (char o=0; o < (char)numObs; o++) {	
		candidate.clear();
		ADDCHRS(candidate, a, o);
		candidate += (*nc)->first;
		test = tests.find(candidate);
		if (test == tests.end()) {
		    // The prefix version doesn't exist yet.
		    pos=tests.size();
		    tests[candidate] = pos;
		    maxTestLen = std::max(maxTestLen, candidate.length());
		}
	    }
	}
    }
    

    // Iterate over test set till it's complete.
    // Add any missing length-2 tests
    // Should not be any observation tests
    // can only add shorter or equal length tests.
    do {
	added=false;
	for (TestsMap::iterator t=tests.begin(); t != tests.end(); t++) {
	    
	    if (t->first.length() == 0) continue;

	    // Check for missing length-2 test
	    candidate = t->first.substr(0, t->first.length() - 2);
	    if (tests.find(candidate) == tests.end()) {
		pos=tests.size();
		tests[candidate] = pos;
		added=true;
	    }
	    
	    // Check for missing observation tests
	    for (char o=0; o < (char)numObs; o++) {
		candidate = t->first.substr(0, t->first.length()-1);
		
		// This one is a suffix test.
		candidate.push_back(o + '0');
		//cout<<"obs cand="<<candidate<<endl;
		if (tests.find(candidate) == tests.end()) {
		    // Add the test
		    pos=tests.size();
		    tests[candidate] = pos;
		    added=true;
		}
	    }
	}
    } while(added);


    testPredictions.resize(maxHistory, tests.size(), true); 
    // True preserves current data

    // init new test predictions to 1.0 then normalise
    project(testPredictions, 
	    range(0, maxHistory), 
	    range(initSize, tests.size())) 
	= scalar_matrix<double>(maxHistory, 
				tests.size() - initSize, 1.0
				);
    
    normalise(row);

}


/**
 * This is the top level routine in the step to step update of the PSR
 * given an incoming observation.
 */
void PSRTransform::getAction(Observation& obs, Vector& action, bool computeGrad) {
    
    // This PSR code can't handle multi-dim obs.
    assert(obs.getFeatures().size1() == 1);
    size_t row = 0;

    // Append most recent history.
    ADDCHRS(history, prevAction, (int)obs.getFeatures()(0, obs.getAgent()));

    // First two parts of history are no longer look ahead history,
    // they are now.  But we might have to add them back in if we grow
    // the maximum test length.  Which will be equivelent to
    // re-updating the same row twice I guess.
    size_t actionDigit =  (size_t)(history[0] - '0');
    size_t obsDigit =  (size_t)(history[1] - '0');
    history.erase(0,2);
    
    // Can't do updates till we've got enough history
    if (history.length() == maxTestLen) {

	// Figure out which row we are about to update. If tests got
	// added, this could go backwards! Some confusion caused
	// because history length is measures as characters, two per
	// step.  obs.steps starts counting at 1, but we already
	// filled in row 0 so that's okay.
	row = ((obs.getSteps()) - maxTestLen/2 + maxHistory)%maxHistory;
	
	// Update history will use 'history' var as a lookahead for 
	// grad ascent
	updateHistory(row, actionDigit, obsDigit);

	// Beware that updating the core tests might cause 
	// problems for innapropriate controllers.
	if (discover && row==maxHistory-1 && row != lastRow && passes++ > passesBeforeAdd) {
	    addCoreTest(row);
	    passes=0;
	}
	
    }


    psrObs.setAgent(obs.getAgent());
    psrObs.setSteps(obs.getSteps());

    // Wind the PSR forward to the current observation
    if (!learn) assert(history.length()==0 && maxTestLen==0);
    windPSRForward(row, history);
    updatePSR((row + history.length()/2)%maxHistory);
    makePSRObs(obsDigit);    
    if (obs.getSteps()%PSR_CONF_FREQ == 0) {
	cout<<"PSRConfidence: "<<confidence<<endl;
	cout<<"PSR: "<<psrVector<<endl;
    }
    try {
	controller->getAction(psrObs, action, computeGrad);
    } catch (exception& e) {
	cout<<"!! "<< e.what()<<endl;
	cout<<"psr: "<<psrVector<<endl;
	cout<<"action: "<<action[0]<<endl;
	cout<<"maxParam: "<<controller->getMaxParam()<<endl;
	abort();
    }

    // The PSR wants to know what action created what observation, not
    // what observation generates which action.
    prevAction = (size_t)action[0];
    assert(prevAction < numActions);
    
    // If maximum test length has grown, we actually need to increase
    // the delay between the current step and the row we are
    // updating. This will mean re-updating the current row, so
    // pre-pend the actionDigit and obsDigit again
    if (learn && history.length() < maxTestLen) {
	string tmpStr;
	ADDCHRS(tmpStr, actionDigit, obsDigit);
	history = tmpStr + history;
    }
    
    lastRow = row;
}


/**
 * Add the previous action and most recent observation to the
 * history, put a new row in the test prediction matrix for the new
 * history.
 * @param row of testPrediction matrix to add to
 * @param current observation
 * @param previous action
 * @param Do least squares to update entries we otherwise cannot.
 *        defaults to false
 */
void PSRTransform::updateHistory(size_t row,  size_t action, size_t obs, bool doLS) {


    TestList lst;    // least squares needed for these tests

    testPredictions(row, 0) = 1.0; // Null test always succeeds.
    assert(action >= 0 && action < numActions);
    assert(obs >= 0 && obs < numObs);


    if (learn) {
	lst.clear();
	// Update each column (test) of the matrix
	for (TestsMap::iterator t=tests.begin(); t != tests.end(); t++) {
	    
	    // Is obs,action,test also a test?	
	    if (!computeEntryFromPredictions(row, t, action, obs)) {
		// No, we have to solve a least squares problem in second round
		lst.push_front(t);
	    }
	    
	}
	
	// At this point we are guarnateed an estimate of the current PSR.
	// And to solve least squares problem for other tests we need the
	// current estimate
	computeLeastSquares(row, lst, false); // false=update row directly
    
	// make sure we are normalised before hillClimb
	normalise(row);
	
	// Do gradient steps for eligible entries
	hillClimbPSR(row, history);
    }
    else updateEntriesFromParams(row, action, obs);
    
    // Always normalise to avoid numerical errors aggregating.
    normalise(row);
}

				 
/**
 * Some tests cannot be updated from just entries in the previous
 * row. These tests are updated by least squares, finding m to
 * minimise ||P(Q|H) m - P(t|H)||. We do this for one of two 
 * reasons, to get update a test, or to compute the PSR
 * parameters which, in future allow computation of all
 * core test updates from only the previous core-test row.
 * During learning we always estimate from least squares.
 * After learning we call this once, then just update from
 * there.
 */
void PSRTransform::computeLeastSquares(size_t row, 
				       TestList& tl,
				       bool keepParams) {
    
    int c=0;
    normalise(row); // We've just come out of estimation from
                    // historical data.
    
    // Create core tests prediction matrix just once for all
    // tests to do LS on.
    for (TestsMap::iterator q=coreTests.begin(); q != coreTests.end(); q++) {
	column(Q, c++) = column(testPredictions, q->second);
    }

    if (!keepParams) updatePSR(row); // For dot products
    
    // Loop over tests to be solved for
    for (TestList::iterator t=tl.begin(); t != tl.end(); t++) {
	
	testHistory.assign(column(testPredictions, (*t)->second));
   	
	// Solve for this test.
	try {
	    UBlasExtras::linearLeastSquares(Q, testHistory, params);
	} 
	catch (exception& e) {
	    cout<<'.';
	    // Fill in with previous value until we collect some more data.
	    if (!keepParams) testPredictions(row, (*t)->second) = 
		testPredictions(PREVROW(row), (*t)->second);
	}
	
	if (!keepParams) {
	    
	    // Now the estimate of the prediction is current PSR state vector
	    // inner prod with params.

	    testPredictions(row, (*t)->second) = inner_prod(params, psrVector);
		       
	}
	else {
	    // Just keep the parameters for later
	    column(allParams, (*t)->second) = params;
	}
    }
}


/**
 * We have enough data in the table to update a test prediction
 * directly from table data, new obs and new action.  
 * @param row that we are updating 
 * @param test that we are updating 
 * @param observed symbol 
 * @param action that was then observed.
 */
bool PSRTransform::computeEntryFromPredictions(size_t row, 
					       TestsMap::iterator t, 
					       size_t action, 
					       size_t obs
					 ) {

    TestsMap::iterator numerator;
    TestsMap::iterator denominator;
    string neededTest;
    string base;

    // We've already set empty test=1.0
    if (t->first.length()==0) return true; 

    ADDCHRS(base, action, obs);
    neededTest = base + t->first;

    numerator = tests.find(neededTest);

    if (numerator == tests.end()) {
	// This can't be a core test or something
	// has gone wrong with test set extension
	assert(coreTests.find(t->first) == coreTests.end());
	return false; // Can't compute entry this way
    }

    denominator = tests.find(base);
    assert(denominator != tests.end()); // The obs/action pair should
					// always be a valid test
					// since it's added in
					// constructor.
    
    testPredictions(row, t->second) = 
	testPredictions(PREVROW(row), numerator->second)/
	testPredictions(PREVROW(row), denominator->second);

    return true;
}


/**
 * Should core tests be updated
 */
void PSRTransform::setDiscover(bool discover) {
    this->discover = discover;
}


/**
 * PSRs must maintain correct probability distributions
 * @param The row that has just been modified
 */
void PSRTransform::normalise(size_t row) {

    std::map<string, double> norms; // Holds normalising terms
    std::map<string, double>::iterator n;
    string normKey;
    string parentKey;

    for (TestsMap::iterator t=tests.begin(); t != tests.end(); t++) {

	// 0 values could cause div by 0 error if that prediction is wrong.
	if (testPredictions(row, t->second) <= PG_MACHINE_EPS) {
	    testPredictions(row, t->second) = PG_MACHINE_EPS;
	}

	// Also update the norms where possible.
	
	// Null test has normalisation factor of 1
	if (t->first.length()==0) norms[t->first] = 1.0;
	else {
	    // Otherwise we strip off the final observation to compute
	    // a key for normalisation factor for this test+action
	    // (See Bowling&McCraken) paper if this is confusing)
	    normKey = t->first.substr(0, t->first.length() - 1);
	    
	    n = norms.find(normKey);
	    if (n == norms.end()) {
		// New entry for this norm
		norms[normKey] = testPredictions(row, t->second);
	    }
	    else {
		// Add to the norm sum.
		norms[normKey] += testPredictions(row, t->second);
	    }
	}
    }


    // Now, in a single pass over the data attempt to normalise all
    // tests in groups that only differ by final observation.
    for (TestsMap::iterator t=tests.begin(); t != tests.end(); t++) {
	
	// 0 length case does not need normalisation because that's 
	// the null test.
    
	if (t->first.length() == 0) continue;
	
	normKey = t->first.substr(0, t->first.length() - 1);
	parentKey = t->first.substr(0, t->first.length() - 2);
	testPredictions(row, t->second) *= testPredictions(row, tests[parentKey])/norms[normKey];
    }   
}


/**
 * Copy the relevant row of the test prediction matrix into 
 * the PSR feature vector. Only copy for tests that are part of the
 * core test set
 * @param The current row of the history matrix.
 */
void PSRTransform::updatePSR(size_t rowv) {

    size_t qIndex=0;
  
    for (TestsMap::iterator q=coreTests.begin(); q != coreTests.end(); q++) {
	if ( isnan(testPredictions(rowv, q->second)) || 
		   testPredictions(rowv, q->second) < 0.0 ||
		   testPredictions(rowv, q->second) > 1.0) {
	    cout<<"row="
		<<rowv
		<<" out of bounds: "
		<<row(testPredictions, rowv)
		<<endl;
	}
	psrVector(qIndex++) = testPredictions(rowv, q->second);
    }
}


/**
 * Potentially move a current test to the status of core test.
 * Performed by singular value decomposition followed by taking the
 * current most likely non-core test.
 */
void PSRTransform::addCoreTest(size_t row) {
    
    assert(row == maxHistory-1);    
    
    if (coreTests.size() == maxCores) {
	cout<<"Core tests full. Switch to PSR params\n";
	
	estimateAllParams();
	discover=false;
	learn=false;
	
	// Catch up to the current observation.
	windPSRForward(row, history);
	size_t newRow = (row + history.length()/2)%maxHistory;
	updatePSR(newRow);
	//cout<<"PSR after final wind forward to row="<<newRow<<" :"<<psrVector<<endl;
	maxTestLen=0;
	history.clear();
	
	// Turn on controller learning
	controller->setStepSize(algStepSize);

	return;
    }
    

    double minCondNum = minKappaForCore;
    double condNum;
    TestList addList;
    int q;

    Matrix corePlus1(maxHistory, coreTests.size() + 1);
    Vector singularValues(coreTests.size() + 1);

    if (tests.size() > maxHistory) {
	throw std::runtime_error("SVD requires at least as much history as tests\n");
    }

    q=0;
    // Load the standard core tests history into the matrix to do SVD on.
    for (TestsMap::iterator c = coreTests.begin(); c != coreTests.end(); c++) {
	column(corePlus1, q++) = column(testPredictions, c->second);
    }


    // Add the next most significant test not already in core set.
    for (TestsMap::iterator t=tests.begin(); t != tests.end(); t++) {

	q = coreTests.size();

	// Ignore current core tests.
	if (coreTests.find(t->first) != coreTests.end()) continue;
 
	// Ignore tests that are not one step prefix extensions
	// i.e., must be in set X.
	if (coreTests.find(t->first.substr(2, t->first.length()-2)) == coreTests.end()) 
 	    continue;

	column(corePlus1, q) = column(testPredictions, t->second);
	UBlasExtras::svd(corePlus1, singularValues);
	
	while (singularValues[q] <= 0.0) q--;
	assert(q>0);
	condNum = singularValues[0]/singularValues[q];
	cout<<"Test '"
	    <<t->first
	    <<"' has kappa="
	    <<condNum
	    <<" svs:"
	    <<singularValues;

	if (condNum < minKappaForCore && 
	    addAllCores && 
	    (coreTests.size()+addList.size()) <  (size_t)controller->getInputDim() ) {
	    addList.push_back(t);
	}    
	else if (condNum < minCondNum) {
	    minCondNum = condNum;
	    addList.clear();
	    addList.push_back(t);
	    //cout<<"new minimum="<<minCondNum<<" '"<<t->first<<"'\n";
	}
	cout<<endl;
    }

    // Add all the new core tests.
    if (!addList.empty()) {
	for (TestList::iterator toa=addList.begin(); toa != addList.end(); toa++) {
	    
	    cout<<"Adding `"<<(*toa)->first<<"'"<<endl;
	    
	    // Actually add the test.
	    coreTests[(*toa)->first]=(*toa)->second;
	   
	    // Add new tests that are proposed extensions of this one.

	}
	completeTestSet(row, addList);
	params.resize(coreTests.size());
	Q.resize(maxHistory, coreTests.size());
	psrVector.resize(coreTests.size());
    }
    else cout<<"No core candidates found :(\n";

}


/**
 * A simple TD based hill climb to reduce the PSR vector error. I.e.,
 * learning only needs to know the row & current history vector to see
 * if tests will come true.  Retains normalisation property
 * @param which row are we operating on.
*/
void PSRTransform::hillClimbPSR(size_t row, string hist) {

    TestsMap::iterator parentCol;
    TestsMap::iterator normTest;
    TestList normList;
    double parentVal;
    string parentTest;
    string normKey;
    string tmpKey;
    double normVal;
    double oldPrediction;
    double desiredVal;

    confidence *= ACCURACY_DISCOUNT;
    // Go though all the tests and search for matches with the suffix
    // of the current history.
    for (TestsMap::iterator t=tests.begin(); t != tests.end(); t++) {

	// not interested in first test.
	if (t->first.length() == 0) continue;

	if (t->first == hist.substr(0, t->first.length())) {
	    // We have a match. Well done little PSR
	    // Now get the parent test prediction
	    
	    parentTest = t->first.substr(0, t->first.length()-2);
	    parentCol = tests.find(parentTest);
	    assert(parentCol != tests.end());
	    parentVal=testPredictions(row, parentCol->second);
	    
	    // Update the prediction.
	    oldPrediction = testPredictions(row, t->second);
	    desiredVal = (1.0 - stepSize)*oldPrediction + stepSize*parentVal;
	    testPredictions(row, t->second) = desiredVal;

	    // Update the crappy but cheap confidence metric
	    confidence += (1 - ACCURACY_DISCOUNT)*oldPrediction;

	    // Now renormalise other tests that are the same except for 
	    // the last observation. This is necessary to get the
	    // correct value for longer tests that might be about
	    // to be updated in hill climb, and can be done
	    // more efficiently here than calling normalise repeatedly.
	    normList.clear();
	    normVal=(parentVal - desiredVal)/
		computeNorm(row, t->first, normList);
	   
	    for (TestList::iterator normTerm = normList.begin(); 
		 normTerm != normList.end(); 
		 normTerm++
		 ) {
		testPredictions(row, (*normTerm)->second) *= normVal;
	    }    
	}
    }
}


/**
 * A testing routine only. Computes likelihood of all predicted next
 * symbols and compares with actual from history.
 */
void PSRTransform::updatePrediction(bool print) {

    map<char, double> likelihoods;
    double maxVote = 0;
    char votedSymbol = 'x';

    if (print) cout<<" Hist='"<<history<<"' PSR="<<psrVector<<endl;
    
    for (TestsMap::iterator q=coreTests.begin();  q != coreTests.end(); q++) {
	// Only take one step tests where the action matches.
	if (q->first.length()==2 && q->first[0] == history[0]) {
	    likelihoods[q->first[1]] += psrVector(q->second); 
	}
    }
    
    if (print) cout<<"Actual="<<history[1]<<" Predictions: ";
    for (map<char, double>::iterator o=likelihoods.begin(); o != likelihoods.end(); o++) {
	if (print) cout<<" '"<<o->first<<"'="<<o->second;
	if (o->second > maxVote) {
	    maxVote = o->second;
	    votedSymbol = o->first;
	}
    }

    predictionAccuracy *= ACCURACY_DISCOUNT;
    if (votedSymbol==history[1] && maxVote>0.0) {
	predictionAccuracy += 1.0 - ACCURACY_DISCOUNT;
    }
    if (print) cout<<" Accuracy="<<predictionAccuracy<<endl;

}


/**
 * Compute the sum of all predictions for all tests 
 * that differ from the parameter by at most the last 
 * observation. Used for normalising predictions.
 * during the hill climb.
 * @param row to operate on
 * @param test to compute normalising terms on.
 * @param empty test list which will be loaded with tests
 * contributed to the norm. Does not include the test itself
 * @return the sum of predictions for all tests that differ
 * by only last obs.
 */
double PSRTransform::computeNorm(size_t row, const string& test, TestList& tl) {

    string normKey = test.substr(0, test.length()-1);
    string srchKey;
    double norm=0.0;
    TestsMap::iterator foundTest;

    for (char o='0'; o < (char)('0' + numObs); o++) {
	srchKey = normKey + o;
	// Don't include the test itself.
	if (srchKey != test) {
	    foundTest = tests.find(srchKey);
	    assert(foundTest != tests.end());
	    norm += testPredictions(row, foundTest->second);
	    tl.push_back(foundTest);
	}

    }

    return norm;
}


/**
 * Wind the clock forward to the current observation, updating the PSR
 * vector appropriately. We can't afford to update all the tests that
 * don't have an extended version in the set, i.e., tests that we have
 * to do regression on. So we ignore those, copying their previous
 * values.
 * @param current history
 * @param current row
 */
void PSRTransform::windPSRForward(size_t row, string hist) {

    size_t actionDigit;
    size_t obsDigit;

    //cout<<"Max test len="<<maxTestLen<<" row="<<row<<" winding forward '"<<hist<<"' current PSR:"<<psrVector<<endl;
    while (hist.length() > 0) {

	// First two parts of history are no longer look ahead history,
	// they are now.  But we might have to add them back in if we grow
	// the maximum test length.  Which will be equivelent to
	// re-updating the same row twice I guess.
	actionDigit =  (size_t)(hist[0] - '0');
	obsDigit =  (size_t)(hist[1] - '0');
	hist.erase(0,2);

	// Get to next row
	row = (row + 1)%maxHistory;
	
	// Update history, but avoid expensive Least squares calcs.
	updateHistory(row, actionDigit, obsDigit, false);

	updatePSR(row);
	//cout<<"PSR winding forward to "<<row<<" :"<<psrVector<<endl;
    }

    // Update getFeatures() with PSR from current row 

    //cout<<"Done wind forward to row="<<row<<" current PSR:"<<psrVector<<endl;

}



void PSRTransform::estimateAllParams() {

    allParams.resize(coreTests.size(), tests.size());
    allParams.clear();

    TestList relevantTests;
    string extendedTests;

    // Relavant tests are all core tests plus one-step prefix
    for (TestsMap::iterator ct=coreTests.begin(); ct != coreTests.end(); ct++) {

	cout<<"estimating params for '"<<ct->first<<"'\n";

	relevantTests.push_back(ct);

	// Find all the one-step prefixe tests. Need their params too.
	for (char a=0; a < (char)numActions; a++) {
	    for (char o=0; o < (char)numObs; o++) {
		extendedTests.clear();
		ADDCHRS(extendedTests, a, o);
		extendedTests += ct->first;
		TestsMap::iterator et = tests.find(extendedTests);
		assert(et != tests.end());
		relevantTests.push_back(et);
	    }
	}
    }
    
    // 0 is dummy row.. how ugly, true=keep params
    computeLeastSquares(0, relevantTests, true); 

    cout<<"Finished estmation: "<<allParams<<endl;

    // Okay, from now on we should be able to update the PSR vector
    // purely from the parameters and previous row.
}


/**
 * Update core test entries for current row based on previous core test values
 * and parameters previously estimates.
 * @param Assumes allParams is up to date
 * @param Current row (will figure out previous row)
 * @param next action on history stack
 * @param next obs on history stack.
*/
void PSRTransform::updateEntriesFromParams(size_t urow, size_t action, size_t obs) {

    //cout<<"Updating row="<<urow<<" from params"<<endl
    //<<"Prev row="<<row(testPredictions, PREVROW(urow))<<endl
    //<<"action="<<action<<" obs="<<obs<<endl;

    // Get previous row
    updatePSR(PREVROW(urow));
    
    TestsMap::iterator neededCoreTest;
    string base;
    string neededCoreString;
    ADDCHRS(base, action, obs);
    size_t baseColumn = tests[base];

    double numerator;
    double denominator;

    for (TestsMap::iterator ct=coreTests.begin(); ct != coreTests.end(); ct++) {
    
	// Compute test with correct one step prefix for this core test
	neededCoreString = base + ct->first;
	neededCoreTest = tests.find(neededCoreString);
	assert(neededCoreTest != tests.end());
	
	// At this point we must have the full core test compliment.
	assert(psrVector.size() == coreTests.size());

	// sanity check that params have been filled in
	assert(norm_inf(column(allParams, neededCoreTest->second)) > 0.0);
	assert(norm_inf(column(allParams, baseColumn)) > 0.0);

	// Remember that psrVector now holds previous values.
	numerator = inner_prod(psrVector, 
			       column(allParams, neededCoreTest->second));
	denominator = inner_prod(psrVector, 
				 column(allParams, baseColumn));

	// Update the entry.
	testPredictions(urow, ct->second) = numerator/denominator;

    }
    
    
}


void PSRTransform::setStepSize(double stepSize) {
 
    // Overide step size until PSR learning has finished.
    algStepSize = stepSize;
    controller->setStepSize(0.0);

}


/**
 * Turn the PSR into a binary format by thresholding on 0.5
 * Want to see if this makes learning easier. 
 */
void PSRTransform::makePSRObs(int obs) {
    
    psrObs.getFeatures().clear();
    for (size_t i=0; i < psrVector.size(); i++) {
	//cout<<"Test
	if (binaryObs) {
	    psrObs.getFeatures()(i, 0)=(psrVector(i)>=PSR_BINARY_THRESH)?1.0:0.0;
	}
	else {
	    psrObs.getFeatures()(i, 0) = psrVector(i);
	}
    }

    if (appendObs) {
	psrObs.getFeatures()(maxCores + obs, 0) = 1.0;
    }

}
}
