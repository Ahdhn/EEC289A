/**
* $Id: GOALPomdp.cc 127 2007-09-10 17:20:04Z daa $
*/


#include<time.h>
#include<iomanip>

#include"PGBasics.hh"

#include"OLPomdp.hh"
#include"GOALPomdp.hh"


using namespace std;

namespace libpg {

	/**
	* Top level piece of glue for Policy gradient, link the controller
	* and the simulator
	* @param  Controller
	* @param  Simulator (problem specific)
	* @param  discount factor [0, 1]. If 1, trace must be reset
	* occasionally or have 0 mean reward.
	*/
	GOALPomdp::GOALPomdp(Controller* controller, Simulator* simulator, double discount, double stepSize) : 
		OLPomdp(controller, simulator, discount, stepSize) {
		episodeSteps = 0;
		totalEpisodes = 0;
	}



	/**
	* Do a fixed number of individual steps, aggregating the reward.  No
	* output is done. Could also be termed an "epoch"
	* @param  empty vector to accumulate average reward for each reward dimension
	* @param  number ofsteps to iterate for
	* @param  should learning be turned on, or are we
	* just exploiting/evaluating the current policy?
	*/
	double GOALPomdp::doSteps(Vector& totalRewards, int steps, bool learn) {

		Vector action(simulator->getActionDim());
		Vector rewards(simulator->getRewardDim());
		Vector reinforcement(simulator->getRewardDim());
		bool goalReached = false;

		totalRewards.clear();
		episodes = 0;
		successes = 0;

		for (int s = 0; s < steps; s++) {

			simulator->getObservation(obs);
			if (learn && obs.getSteps()>0 && !goalReached) {
				// Must do next line before calling getAction() since
				// getAction will accumulate log action gradients into the
				// trace directly for efficiency.
				controller->discountTrace();
			}
			goalReached = false;
			controller->getAction(obs, action, learn);
			episodeSteps++;
			obs.setSteps(obs.getSteps() + 1);
			goalReached = (bool)simulator->doAction(action);
			simulator->getReward(rewards);
			totalRewards += rewards;
			if (learn) {
				reinforcement.assign(rewards - baseline);
				// [daa] Check for 0 rewards is done at Controller level
				controller->instantStep(reinforcement);
			}
			if (goalReached) {
				controller->resetTrace();
				episodeSteps = 0;
				episodes++;
				// A positive reward counts as a successful run
				if (norm_1(rewards) > 0) successes++;
				rewards.clear();
			}
		}

		totalEpisodes += episodes;
		totalRewards /= steps;
		return inner_prod(totalRewards, scalar_vector<double>(simulator->getRewardDim(), 1.0)) / simulator->getRewardDim();
	}


	/**
	* Do as many steps as required to get to the goal for the first time.
	* Assumes simulator has been reset. EpisodeSteps will be reset to 0.
	* @param empty vector to accumulate average reward for each reward dimension
	* @param should learning be turned on, or are we
	* just exploiting/evaluating the current policy?
	* @returns average reward. Can also query class for episodeSteps
	*/
	double GOALPomdp::doSingleEpisode(Vector& totalRewards, bool learn) {

		Vector action(simulator->getActionDim());
		Vector rewards(simulator->getRewardDim());
		Vector reinforcement(simulator->getRewardDim());
		bool goalReached = false;

		episodeSteps = 0;
		successes = 0;
		episodes = 0;

		totalRewards.clear();

		do {

			simulator->getObservation(obs);
			if (learn) {
				// Must do next line before calling getAction() since
				// getAction will accumulate log action gradients into the
				// trace directly for efficiency.
				controller->discountTrace();
			}
			controller->getAction(obs, action, learn);
			episodeSteps++;
			obs.setSteps(obs.getSteps() + 1);
			goalReached = (bool)simulator->doAction(action);
			simulator->getReward(rewards);
			totalRewards += rewards;
			if (learn) {
				reinforcement.assign(rewards - baseline);
				// [daa] Check for 0 rewards is done at Controller level
				controller->instantStep(reinforcement);
			}
		} while (!goalReached);

		controller->resetTrace();

		totalRewards /= episodeSteps;

		episodes++;
		if (norm_1(rewards) > 0) successes++; // A positive reward counts as a successful run

		// Return a scalar performance measure from vector rewards.
		return inner_prod(totalRewards, scalar_vector<double>(simulator->getRewardDim(), 1.0)) / simulator->getRewardDim();
	}


	/**
	* Print a line of performance information
	* @param  are we only printing the headers?
	* @param  total number of steps so far
	* @param  average reward over last epoch
	* @param  maximum parameter
	* @param  time elapsed so far.
	* This overloaded method from RLAlg additionally prints prob(success)
	*/
	void GOALPomdp::printPerformanceInfo(bool printTitles, int steps, double avgReward, double maxParam, int seconds) {
		cout << setprecision(PRECISION);
		if (printTitles) {
			cout << "       "
				<< setw(CWIDTH) << "Steps"
				<< setw(CWIDTH) << "R"
				<< setw(CWIDTH) << "max(Param)"
				<< setw(CWIDTH) << "Seconds"
				<< setw(CWIDTH) << "Episodes"
				<< setw(CWIDTH) << "Succ%" << endl;
		}
		cout << "Trial: "
			<< setw(CWIDTH) << steps
			<< setw(CWIDTH) << avgReward
			<< setw(CWIDTH) << maxParam
			<< setw(CWIDTH) << seconds
			<< setw(CWIDTH) << totalEpisodes
			<< setw(CWIDTH) << successes / (double)episodes << endl;
	}


}
