Policy-Gradient Library
=======================

Revision 
This code is distributed under the Mozilla Public License V1.1. 
See MPL-1.1.txt

Contents:

1. Overview
2. Main classes
3. Installation
4. Simulators
5. Getting Started
6. Contributors
7. Known Flaws
8. Wishlist

******************************************************************** 
IMPORTANT NOTES:

The full API is documented in docs/html/index.html 
If you have checked the source out of a repsository, you need to run
doxygen docs/Doxyfile to generate the API docs.

Worked example: The ThreeState simulator is a really simple MDP that
demonstrates most of the algorithms in libpg. It is in the 
simulators directory, but there is an accompanying document in
docs/WorkedExample.

The quickest way to code up your optimisation problem for solving
with RL is to fill in the blanks in the template code in 
libpg/simulators/TemplateRLApp. 
********************************************************************


1. Overview
-----------

The PG library was intended to be a high-performance policy-gradient
reinforcement learning library. Since the first version it has been
extended to a number of value based RL algorithms, so the name is only
historical. It's now a general RL library.  It has been designed with
large distributed RL systems in mind. It's not perfect, but it is
pretty fast. This README contains all the instructions you are going
to get on building and compiling, but see the Doxygen files for
overviews of each class and their APIs. The ThreeState example
simulator has also been constructed to show how to invoke almost all
the RL algorithms in this library on a simple 3 state MDP simulator.

What libpg does NOT provide is model based planning algorithms such as
value iteration, or real-time dynamic programming, or exact policy
gradient. There is limited support for belief state tracking in the
simulators/Cassandra/ directory (named because we use the POMDP file
format created by Anthony Cassandra).  One day I'd like to extend it
to these situations, but that will require some uptake of the library.

1.1 Project goals:

o Provide easy to use implementations of state-of-the-art RL algorithms 
  for the non RL savvy, allowing immediate application to difficult industry
  problems. 
o High performance, especially on multi-agent RL problems
o Extensible plug'n'play algorithms for research purposes

1.2 Main algorithms implemented:

RL algorithms:
o SARSA
o QLearning
o Vanilla Policy-Gradient (online or batch, including line search)
o Natural actor-critic 
o Least Squares TD-Q(\lambda) (for LSPI)
o Least Squares Policy Iteration

Misc supporting algorithms:
o Line searches for batch mode
o Tikhonov Regularisation
o HMM based POMDP state estimation
o Finite history transformation
o Importance sampling and bias for including hard wired policies

In specific simulators:
o Predictive state representation learning and transformations
o Belief state tracking
o Conditional random fields for DEC-POMDPs (i.e., efficient Collaborative RL)
o Parser for Anthony Cassandra's POMDPs.

There are six important abstract classes:

2.1 Approximator:

An approximator provides a parameterised function approximator,
mapping observations to action likelihoods (not well formed
probability distributions). For speed and ease of creating distributed
RL systems where each agent can have its own parameters/discount
factors/rewards/whatever, much of the algorithm implementation is
pushed all the way down to the approximators. I.e., approximators know
how to make gradient steps and accumulate eligibility traces.

  Examples:
	LookupTable: provides tabular rasa mapping
	LookupTableBatch: extension to work with batch algorithms
	NeuralNet: provides a simple neural network
	NeuralNetBatch: extension to batch algs
	NeuralNetAtlas: uses the ATLAS library for fast neural nets.

2.2 Controller:

Controllers serve multiple purposes, and complex distributed RL
applications can be created by chaining together sequences of
Controllers.

A controller provides the interface between the learning algorithm and
the approximator. Controllers turn the output of Approximators into a
probability distribution, typically through a Gibbs
distribution. Controllers therefore also know how to compute the log
action probability gradient with respect to the outputs of the
Approximator.
  
Unless obvious from the name, classes implement policy-gradient algorithms.

  Examples:
	BasicController: choose action from n possibilities
	BinaryController: optimised for 2 actions only
	ValueController: abstract class for value based methods where 
			 the approximators are assumed to be 
			 approximating values. This class is made
			 concrete by SARSAController, QLearningController,
			 and LSPIController.

Controllers are also the mechanism for achieving distributed RL with
local or non local rewards. 

  Examples:
	FactoredController: groups many controllers, one action per agent.
	ManyToOneController: groups many controllers that individually 
                             vote for their action to be picked over 
                             that of the other controllers.

Controllers can also transform observations to new approximator features:

  Examples:
	NACTransform: implements natural actor-critic PG by linear
		      transforms on features and gradients.
	PSRTransform: implements predictive state representation
	MemoryTransform: implement finite memory window reprensetation

2.3 Simulator:

This implements your domain. It needs to provide Observations,
implement actions, and return rewards. After that, it's up to
you. Basic example is in simulators/ThreeState/ThreeState.cc

2.4 Observation:

A glorified data structure that passes observations
around. Observations can be an arbitrary matrix, but typically each
column in the matrix would be the observation for a particular
agent. The observation also knows how many steps have been taken so far
for algorithms that need to count iterations. It also knows which
agent it is being passed to, so could be used to implement inter-agent
explicit communication if desired.

2.5 RLAlg:

RLAlgs tie everything together, implementing the high level algorithm
by calling appropriate simulator and controller methods. There's
always a single simulator and a single controller at the top level
(although the top level controller can implement a tree of
sub-controllers for distributed apps.)

  Examples: OLPomdp: a generic stochastic gradient ascent. With the
	             right TransformController can be used to
	             implement Baxter Bartletts PG, Natural gradients,
	             Natural actor-critic and PSR learning.
            
	    GPomdp:  Compute batch gradient estimates, but still a 
	    	     fixed step size.
		     
            LineSearchAlg: Batch gradient estimates with a line search
                           designed by Baxter and Bartlett for noisy
                           gradients.

	    GOALPomdp: OLPomdp tailored for shortest stochastic path 
	               problems.

All these algorithms are invoked with the "learn" methods, which is
the same across all agorithms (so far), which gives you print outs of
progress. No other part of the code gives messages unless something
has gone wrong.

The actual work is done by methods doSteps() and
learnCore(). doSteps() does actual iterations, and in online algorithms
does parameter updates. learnCore() calls doSteps for a fixed number
of steps then does batch parameter updates.


2.6 Policy

Values based methods (see ValueController) require an explicit policy 
that maps values to actions. You can subclass Policy to make your own
or use:

  Examples: eGreedyPolicy: epsilon greedy. Acts greedily most of the time, 
  	    		   or uniform randomly the rest of the time.

            SoftMaxPolicy: Maps values to a probability distribution. Higher
	    		   values = higher probability. Uses a temperature
                           parameter to control the peakiness. 

3. Installation
---------------

3.1 Dependencies:

    REQUIRED:
    Boost: Needed for the uBlas C++ linear algebra package.

    OPTIONAL:
    uBlas sandbox: Extra Boost package that provides bindings to the following.
    ATLAS: only if you want to use fast neural nets. Needs Boost Sandbox.
    CLAPACK: if you want to do PSRs or anything else that needs 
    	     singular value decomposition. Needs Boost Sandbox.

3.2 Installation:

Tested so far on Linux, OSX, and Windows with Cygwin/gcc.

3.2.1 Copy Makefile.default to Makefile
3.2.2 Edit Makefile to set appropriate paths and defines for optional packages.
3.2.3 '% make' will give you libpg.a

3.2.4 Go into simulators/ThreeState as an example of an application.
3.2.5 Again, copy Makefile.default to Makefile
3.2.6 Again, edit Makefile to point to any extra libraries and turn on
      options you want.
3.2.7 '% make' will give you  working ./ThreeState program.

The makefile is fairly automatic once you've set your necessary
paths. If you have everything in system default paths you might not
even need to edit the makefile. You can also copy
simulators/ThreeState/Makefile.default to other directories to set up
your own application. Just change the name of the default executable
in the Makefile. It will automatically compile all .cc files in the
directory.


4. Simulators
-------------

If you haven't already figured it out, policy-gradient applications are
created by writing at least two pieces of code (same file or different files)

	* A sub class of Simulator that implements your problem
	* A main routine that follows the following general pattern:

	4.1 Create an instance of your Simulator
	4.2 Create appropriate Approximators
	4.3 Create a Basic, Binary, or value based controller, giving it 
	    the Approximator (and a Policy in the value case).
	4.4 Create other controllers if necessary, to implment more
	    complex learners, or multi-agent applications. Each controller	
	    is fed the more basic ones.
	4.5 Create your RLAlg of choice, giving it both simulator and top
	    level controller.
	4.6 Call RLAlg.learn(steps per epoch, discount factor, step size)
	4.7 Watch learning take place.

	Parameters can be saved and restored via the appropriate
	methods in RLAlg. RLAlg can also automatically save
	whenever a the long-term reward reaches a new peak.

5. Getting Started
------------------

There are many algorithms. See libpg/simulators/ThreeStateTest.cc for
examples of how to invoke the main ones.

If you just want to get your RL application running as fast and easily
as possible, fill in the template code in
libpg/simulators/TemplateRLApp/. This will assume you want to use LSPI
or Natural Actor Critic. As far as I know, these are the state of the
art in value based and policy based methods respectively.

Due to the use of the uBlas library, you will get a BIG (10x)
performance improvement by compiling with the flag -DNDEBUG. But
beware, this will turn off all assert statements and turn off matrix
bound checks. Make sure your code works first!

The average reward can be subtracted from the immediate
rewards at each step. The average is the average reward gained over
the last "epoch", i.e., last call of doSteps(). Turn this on in RLAlg
with setAutoBaseline(true). I've yet to see a case where the
baseline hurts, but it can produce confusing results. It does pay to
be aware that you need to balance your epoch size by considering
accurate estimates of the baseline with frequency of updating
it. Also, setAutoBaseline=false can be useful in situations where most
rewards are zero, so that a lot of computations can be avoided
(because there is nothing to reinforce).

6. Contributors:
----------------

Main author: Douglas Aberdeen (doug.aberdeen@gmail.com)
ATLAS and SSP code: Olivier Buffet (olivier.buffet@nicta.com.au)
NAC code: Jin Yu (jin.yu@nicta.com.au)
SARSA, Q-Learning, Policy classes: Fabio Pakk Selmi-Dei
Conditional Random fields: Xinhua Zhang
Traffic (SP38): Tony Lopes

Extra contributions are welcome, but please follow coding style that
has been used so far and provide lots of comments in Javadoc format.
Submit contributions to doug.aberdeen@gmail.com

7. Known Flaws
--------------

Most of these issues are due to the fact that historically the library
started out as a policy-gradient library only.

1. Confusing names: OLPomdp is a basic online RL algorithm (value or policy-grad)
   	     	    GPomdp is a basic batch RL algorith (value or policy-grad)
		    Controller classes BasicController and BinaryController are the
		    	       fundamental policy gradient algorithms.

2. Natural actor critic and LSTDQ (and therefore LSPI) assume single agents only.
   Hopefully this will be fixed in the future.

3. I get "Not Yet Implemented" errors. Typically caused by failing to use a Batch
   version of an approximator, or by getting an overloaded prototype wrong. 

4. Q: Do I use the Batch version of the Approximators or non-batch?
   A: It's always safe to use the Batch version 

5. Your documentation sucks! 
   -- Well, give me feedback and I'll fix it.

6. PSR code is hard to understand and use.
   -- Yes.. again, I need feedback to improve it.

7. Some combinations of controllers are incompatible. Yes, currently
the best way to determine this is to understand the algorithms. If you
don't, then try the combination out and see if it works ;)

8. Memory leaks: There are many places where I have made assumptions
   that a constructur will only get called once during runtime, so I
   don't bother to create a nice destructor to deallocate the
   memory. This might affect you if you have a big loop over many
   thousands of optimisation runs.

8. Wishlist
-----------

MPI based concurrency for policy-gradient algorithms.

Standard implementations of benchmark problems like 
mountain car.
