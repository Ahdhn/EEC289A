#ifndef PSRTransform_hh
#define PSRTransform_hh
// $Id: PSRTransform.hh 127 2007-09-10 17:20:04Z daa $
namespace libpg {

/**
 * Take observations and turn them into PSR prediction vector.
 * Optinally learn PSR online.
 */

class PSRTransform : public TransformController {

private:

    std::string nullTest; // Initial core test. Always succeeds.
    Matrix testPredictions; // Rolling matrix with all history elements
    std::string history; // Current history
    size_t maxTestLen; // Current maximum test length
    size_t maxHistory; // Length of history to work with.
    size_t numActions; // Number of action symbols
    size_t numObs; // Number of observation symbols
    Vector params;   // temporary weights for PSR update
    Vector psrVector;   // temporary weights for PSR update
    Vector testHistory; // rhs for least squares
    Matrix allParams;   // permanent weights for PSR update
    Matrix Q; // lhs for least squares
    bool discover; // Are we learning the core tests 
    /**
                    * online, or just tracking psr?
     */
    bool learn; // Are we learning to track the PSR?
    bool addAllCores; // Should we add
    bool binaryObs;
    bool appendObs;
    double minKappaForCore; // Threshold for adding core tests
    int passesBeforeAdd; // How many times to fill testPredictions
			 // before adding core tests.

    double algStepSize; // Step size to pass down to other contorllers

    int passes; // How many times have we filled history matrix.


    /**
     * This trick slightly changes the rules for string ordering.
     * Shorter strings must always preceed longer strings.  This means
     * that the prediction normalisation is done in the correct order
     * just using the normal map iterator.
     */
    struct testEq {
	bool operator()(const std::string& s1, const std::string& s2) const
	{
	    if (s1.length() < s2.length()) return true;
	    if (s1.length() > s2.length()) return false;
	    return s1 < s2;
	}
    };
    
    typedef std::map<std::string, size_t, testEq> TestsMap;

    TestsMap coreTests;
    TestsMap tests; // Mapping of column number in history to test prefixes.

    /**
     * This will prove useful in keeping lists of tests.
     */
    typedef std::list<TestsMap::iterator> TestList;
    
    size_t prevAction; // The action that generated the current obs
    size_t lastRow;
    size_t maxCores;

    Observation psrObs;

    double stepSize; // PSRs have their own gradient step :(

    double predictionAccuracy;
    double confidence;

    void updatePrediction(bool print);
    void updateHistory(size_t row, size_t action, size_t obs, bool doLS=true);
    void updatePSR(size_t row);
    /**
     * Do SVD decomp of History to find new core tests.
     */
    void addCoreTest(size_t row); 

    bool computeEntryFromPredictions(size_t row, 
				     TestsMap::iterator t, 
				     size_t action, 
				     size_t obs);

    void computeLeastSquares(size_t row, 
			     TestList& tl, 
			     bool keepParams);

    void normalise(size_t row);
    
    void completeTestSet(size_t row, TestList& newCores);

    void hillClimbPSR(size_t row, std::string hist);

    double computeNorm(size_t row, const std::string& test, TestList& tl);


    void windPSRForward(size_t row, std::string hist);

    void estimateAllParams();
    void updateEntriesFromParams(size_t row, size_t action, size_t obs);
    void makePSRObs(int obs);

public:
    
    PSRTransform(Controller* controller, 
		 size_t maxHistory,
		 size_t numObs,
		 size_t numActions,
		 double stepSize,
		 int passesBeforeAdd,
		 double minKappaForCore);

    ~PSRTransform() {};

    void setDiscover(bool discover);

    virtual void setStepSize(double stepSize);

    virtual void getAction(Observation& obs, 
			   Vector& action, 
			   bool computeGrad = true);
};
}
#endif
