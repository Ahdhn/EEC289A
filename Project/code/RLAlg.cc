/**
* $Id: RLAlg.cc 127 2007-09-10 17:20:04Z daa $
*/

#include"PGBasics.hh"
#include<time.h>
#include<iomanip>

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
	RLAlg::RLAlg(Controller* controller, Simulator* simulator, double discount, double stepSize) {

		this->controller = controller;
		this->simulator = simulator;
		this->stepSize = stepSize;

		useAutoBaseline = false;

		saveName = NULL;
		doSaves = false;
		epochsBetweenStepUpdates = 0;

		controller->setDiscount(discount);
		controller->setStepSize(stepSize);

		assert(simulator->getRewardDim() != 0);
		baseline.resize(simulator->getRewardDim());
		resetLearning();

		// Initialise the observation class for all possible learners.
		// They can reinitialise if they like, but this should always be right I think.
		obs.init(simulator->getObsRows(), simulator->getObsCols(), simulator->getAgents(), simulator->getNumActions());

	}


	void RLAlg::resetLearning() {

		baseline.clear();
		controller->resetTrace();
		controller->resetParams();
		obs.setSteps(0);

	}


	/**
	* Print a line of performance information
	* @param  are we only printing the headers?
	* @param  total number of steps so far
	* @param  average reward over last epoch
	* @param  maximum parameter
	* @param  time elapsed so far.
	*/
	void RLAlg::printPerformanceInfo(bool printTitles, int steps, double avgReward, double maxParam, int seconds) {
		cout << setprecision(PRECISION);
		if (printTitles) {
			cout << "       "
				<< setw(CWIDTH) << "Steps"
				<< setw(CWIDTH) << "R"
				<< setw(CWIDTH) << "max(Param)"
				<< setw(CWIDTH) << "Seconds" << endl;
		}
		cout << "Trial: "
			<< setw(CWIDTH) << steps
			<< setw(CWIDTH) << avgReward
			<< setw(CWIDTH) << maxParam
			<< setw(CWIDTH) << seconds << endl;
	}


	/**
	* Do optimisation until either time expires or steps runs out. We
	* train in epochs of a fixed number of steps, reporting after each
	* epoch. We get the initial reward by running for an epoch without
	* learning.  This initial epoch does not count towards step totals or
	* final long-term average reward or time
	* @param  steps per epoch
	* @param  time limit (0 for no limit [default])
	* @param  steps (0 for no limit [default])
	* @param  max value (0 for no limit [default])
	* @return  final epoch reward
	*/
	double RLAlg::learn(int stepsPerEpoch, int maxTime, int ms, double maxVal) {
		
		this->maxSteps = ms;
		this->maxTime = maxTime;
		this->maxVal = maxVal;

		int totalSteps = 0;
		int secs = 0;
		int epoch = 0;
		double lastEpochVal;
		double lastStepCheckVal;
		double bestEpochVal;
		Vector totalRewards(simulator->getRewardDim());

		totalRewards.clear();

		time_t startTime = time(NULL);

		bestEpochVal = lastEpochVal = doSteps(totalRewards, stepsPerEpoch, false);
		lastStepCheckVal = bestEpochVal;
		// Record time so far.
		secs = time(NULL) - startTime;
		printPerformanceInfo(true, 0, lastEpochVal, controller->getMaxParam(), secs);

		// Just do this once at the beginning.
		controller->resetTrace();

		while (true) {

			// Halting conditions.
			if (testForHalt(secs, totalSteps, lastEpochVal)) break;

			// Do one learning iteration over stepsPerEpoch
			lastEpochVal = learnCore(stepsPerEpoch, totalSteps);


			// Record time so far.
			secs = time(NULL) - startTime;

			// Save params if better than before
			if (lastEpochVal > bestEpochVal) {
				bestEpochVal = lastEpochVal;
				if (doSaves) write(saveName);
			}

			printPerformanceInfo(false,
				totalSteps,
				lastEpochVal,
				controller->getMaxParam(),
				secs);

			if (epochsBetweenStepUpdates > 0 && epoch%epochsBetweenStepUpdates == 0) {
				checkStepSize(lastEpochVal, bestEpochVal, lastStepCheckVal);
				lastStepCheckVal = lastEpochVal;
			}

			epoch++;
		}

		return lastEpochVal;

	}


	/**
	* Check all the halting conditions. This is provided for convenience of
	* users to override for their own halting conditions.
	* @param wall clock time how long the current learn() call has been running
	* @param steps the total learning steps in the current learn() call
	* @param val the current long term average reward (over last epoch)
	* @return true if optimisation should halt
	*/
	bool RLAlg::testForHalt(time_t secs, int steps, double val) {

		// Halting conditions
		if (maxTime && secs >= maxTime) return true;
		if (maxSteps && steps >= maxSteps) return true;
		if (maxVal && val >= maxVal) return true;
		return false;
	}



	/**
	* By default this is very simple, call doSteps for stepsPerEpoch times,
	* leaving online learning to doSteps. In other RLAlgs this might
	* call doSteps to build up a batch estimate and then do a line search.
	*/
	double RLAlg::learnCore(int stepsPerEpoch, int& totalSteps) {

		double lastEpochVal;
		Vector totalRewards(simulator->getRewardDim());
		totalRewards.clear();

		// Learn
		lastEpochVal = doSteps(totalRewards, stepsPerEpoch, true);
		totalSteps += stepsPerEpoch;
		if (useAutoBaseline) {
			baseline.assign(totalRewards);
		}

		return lastEpochVal;

	}


	/**
	* Provide the filename for regular paramter saves
	* whenever we achieve a new best policy
	*/
	void RLAlg::saveBest(char* fname) {

		saveName = strdup(fname);
		doSaves = true;

	}


	/**
	* Create an output stream and write all controller parameters to it.
	* @param  file name to write to
	*/
	void RLAlg::write(char* fname) {

		ofstream o(fname);
		if (o.is_open()) {
			controller->write(o);
			o.close();
		}
		else cout << "!! Failed to open " << fname << " for write\n";

	}


	/**
	* Restore controller parameters. Be careful to make sure your reading
	* params from the same controller arch.  There's no checks for this
	* at present.
	* @param  file name.
	*/
	void RLAlg::read(char* fname) {

		ifstream o(fname);
		if (o.is_open()) {
			controller->read(o);
			o.close();
		}
		else {
			cout << "!! Failed to open " << fname << " for read\n";
			exit(EXIT_FAILURE);
		}
	}


	/**
	* Do automatic scalar step size update. Will be redundant if using
	* SMD.  Very dumb: attempts to back off quickly if policy dies
	* catastrophically, and slowly increase if we stop learning fast
	* enough. Should help to go to a more deterministic policy.
	* This is not a substitute for a good initial step size, but it
	* will hopefully save you if things go bad half way through optimisation.
	* @param  current average reward
	* @param  best average reward
	* @param  average reward **last time checkStepSize** was called.
	*/
	void RLAlg::checkStepSize(double curr, double best, double last) {

		// We should only be doing this if user requested auto
		// step size updates
		assert(epochsBetweenStepUpdates > 0);

		if (curr < best / (1.0 + threshold)) {
			// We've gotten substantially worse. Panic!
			stepSize /= decreaseFactor;
			controller->setStepSize(stepSize);
			cout << "Policy getting worse. Step size decreased to " << stepSize;
			if (saveName != NULL) {
				// Load the last best parameters if they are around
				read(saveName);
				cout << " and best parameters restored";
			}
			cout << endl;
			return;
		}


		if (curr < last*(1.0 + threshold)) {
			// We're learning too slowly. Speed up.
			stepSize *= increaseFactor;
			controller->setStepSize(stepSize);
			cout << "Policy improving too slowly. Step size increased to " << stepSize << endl;
			return;
		}


	}

	/**
	* Restore the best parameters recorded in file "savename".
	* Should be used after optimisation, before using the policy.
	*/
	void RLAlg::restoreBestParameters(){
		if (saveName != NULL) {
			read(saveName);
			cout << "Best parameters restored";
		}
		else
			cout << "No best parameters to restore !";
		cout << endl;
		return;
	}

	void RLAlg::setStepSize(double s){
		stepSize = s;
	}

	double RLAlg::getStepSize() {
		return stepSize;
	}


	/**
	* Just do an evaluation phase with the most probable action. This
	* gives a different answer from doSteps without learning because
	* doSteps is still following the stochastic policy.  If you are using
	* your own controller, and you haven't implemented getMostProbAction,
	* you will get the default behavious which is to just call getAction.
	* @param  empty rewards vector to populate
	* @param  number of steps to evaluate for
	*/
	double RLAlg::evaluate(Vector&totalRewards, int steps) {

		Observation obs(simulator->getObsRows(),
			simulator->getObsCols(),
			simulator->getAgents(),
			simulator->getNumActions());
		Vector action(simulator->getActionDim());
		Vector rewards(simulator->getRewardDim());

		totalRewards.clear();

		for (int s = 0; s < steps; s++) {
			simulator->getObservation(obs);
			// Must do next line before calling getAction() since
			// getAction will accumulate log action gradients into the
			// trace directly for efficiency.
			controller->getMostProbAction(obs, action);
			simulator->doAction(action);
			simulator->getReward(rewards);
			totalRewards += rewards;
		}

		totalRewards /= steps;

		return norm_1(totalRewards) / simulator->getRewardDim();

	}


	/**
	* Evaluate for a number of steps, then print the results
	* @param  number of steps to evaluate before.
	*/
	double RLAlg::evaluateAndPrint(int stepsPerEpoch, int maxSteps) {
		
		double val;
		Vector totalRewards(simulator->getRewardDim());

		totalRewards.clear();

		time_t startTime = time(NULL);

		//cout<<setprecision(PRECISION);

		val = evaluate(totalRewards, stepsPerEpoch);
		//printPerformanceInfo(true, 0, 0, 0, 0);  
		printPerformanceInfo(true, 0, val, controller->getMaxParam(), time(NULL) - startTime);
		int totalSteps = stepsPerEpoch;
		while (!(maxSteps && totalSteps >= maxSteps)) {
			val = evaluate(totalRewards, stepsPerEpoch);
			printPerformanceInfo(false, totalSteps, val, controller->getMaxParam(), time(NULL) - startTime);
			totalSteps += stepsPerEpoch;
		}
		return val;
	}
}
