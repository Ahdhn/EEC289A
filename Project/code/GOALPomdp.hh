#ifndef GOALPomdp_hh
#define GOALPomdp_hh
// $Id: GOALPomdp.hh 127 2007-09-10 17:20:04Z daa $

namespace libpg {
    class GOALPomdp : public OLPomdp {


    protected:

	int episodes;   // How many episodes in last call to doSteps
	int successes; // How many of those episodes were successful

	int totalEpisodes; // Total episodes over optimisation

	/**
	 * Don't call this one. Tears will result.
	 */
	GOALPomdp() {}; 
    
    public:

	/**
	 * Shut up compiler warning
	 */
	virtual ~GOALPomdp() {};

	int episodeSteps;
    

	GOALPomdp(Controller* controller, Simulator* simulator, double discount, double stepSize);
	virtual double doSteps(Vector& totalRewards, int steps, bool learn);
	double doSingleEpisode(Vector& totalRewards, bool learn);
	virtual void printPerformanceInfo(bool printTitles, int steps, double avgReward, double maxParam, int seconds);
	
	
    };
}
#endif
