#ifndef RLAlg_hh
#define RLAlg_hh
// $Id: RLAlg.hh 127 2007-09-10 17:20:04Z daa $

namespace libpg {
	class RLAlg {

	protected:

		double maxVal;
		time_t maxTime;

		/**
		* Baseline for reward measure. The reward for each step will have
		* the baseline subtracted from it. Idea is to reduce the variance
		* in gradient estimates by turning the rewards into a zero mean
		* process. Typical baseline is the long-term average reward, but
		* that's not known.
		* DEFAULT: 0.0
		* SEE-ALSO: useAutoBaseline();
		*/
		Vector baseline;

		/**
		* If true, baseLine will always be set to the long-term average
		* reward from the last epoch DEFAULT: false
		*/
		bool useAutoBaseline;


		/**
		* Current step size
		*/
		double stepSize;

		/**
		* Maximum number of steps allowed in total. 0=unlimited.
		*/
		int maxSteps;

		/**
		* Shall we automatically update the scalar step size using a
		* naive algorithm.. Use with caution.  Probably still requires a
		* degree of fine tuning for a class of problems.
		* Any SMD algorithm will ignore this.
		*/
		int epochsBetweenStepUpdates; // 0 for no auto updates
		double threshold; // How much better or worse before taking action?
		double increaseFactor; // How much to increase step size by if slow
		double decreaseFactor; // How much to decrease stepsize if bad

		/**
		* If this is true we save whenever we get a best value so far
		*/
		bool doSaves;
		char* saveName;

		/**
		* The top level controller, possibly with many babies under it
		*/
		Controller* controller;

		/**
		* The simulator that knows about the environment
		*/
		Simulator* simulator;

		/**
		* The observations passed to the simulator and controller
		*/
		Observation obs;

	public:

		RLAlg() { };
		RLAlg(Controller* controller, Simulator* simulator, double discount, double stepSize);
		virtual ~RLAlg() {};

		/**
		* The only abstract method
		*/
		virtual double doSteps(Vector& totalRewards, int steps, bool learn) = 0;

		/**
		* Just do an evaluation phase with the most probable action. This
		* gives a different answer from doSteps without learning because
		* doSteps is still following the stochastic policy.
		*/
		virtual double evaluate(Vector&totalRewards, int steps);

		/**
		* Acts a bit like learn, where evaluate is the equivalent of doSteps
		* Does  sets of stepsPerEpoch evaluations, but includes timing and printing
		* of performance information.
		*/
		virtual double evaluateAndPrint(int stepsPerEpoch, int maxSteps);

		/**
		* Automatically estimate the baseline over an epoch and subtract
		* that baseline from any rewards in the next epoch. Do not use if
		* your rewards are mostly 0 (results in less efficiency in
		* gradient update), but generally a good idea otherwise.
		*/
		virtual void setUseAutoBaseline(bool ubl) { useAutoBaseline = ubl; }

		/**
		* Set the maximum number of steps allowed before learning stops.
		* A value of 0 means unlimited steps.
		*/
		virtual void setMaxSteps(int ms) { maxSteps = ms; }

		/**
		* Terminate when a threshold value is reached.
		* 0 means no threshold (bummer if you want 0 to be the threshold)
		*/
		virtual void setMaxValue(int mv) { maxVal = mv; }

		virtual void resetLearning();

		virtual double learn(int stepsPerEpoch, int maxTime = 0, int maxSteps = 0, double maxValue = 0);
		virtual double learnCore(int stepsPerEpoch, int& totalSteps);

		virtual void saveBest(char* fname);
		virtual char* getSaveName() { return saveName; }

		/**
		* Return a pointer to the current controller
		*/
		virtual Controller* getController() { return controller; }

		virtual void write(char* fname);
		virtual void read(char* fname);
		virtual void printPerformanceInfo(bool printTitles,
			int steps,
			double avgReward,
			double maxParam,
			int seconds);
		virtual void checkStepSize(double curr, double best, double last);
		virtual void restoreBestParameters();

		virtual void setStepSize(double s);
		virtual double getStepSize();

		/**
		* User can override halting conditions with this function.
		* @return true if user should halt.
		*/
		virtual bool testForHalt(time_t secs, int steps, double val);

		/**
		* If non-zero attempt to use a naive automatics step size update.
		* You'll need to tune the threshold, increasFactor, and decreaseFactor
		*/
		void setEpochsBetweenStepUpdates(int e) { this->epochsBetweenStepUpdates = e; }

		/**
		* How much better or worse are tings allowed to get before attempting
		* to change the step size.
		*/
		void setThreshold(double t) { threshold = t; }

		/**
		* How much to increase step size by if slow
		*/
		void setIncreaseFactor(double ifa) { increaseFactor = ifa; }

		/**
		* How much to decrease stepsize if bad
		*/
		void setDecreaseFactor(double df) { decreaseFactor = df; }

	};
}
#endif
