/**
 * $Id: NeuralNetAtlas.cc 127 2007-09-10 17:20:04Z daa $
 */


#if (HAVE_BOOST_SANDBOX && HAVE_ATLAS)

#include<boost/numeric/bindings/traits/ublas_matrix.hpp>
#include<boost/numeric/bindings/traits/ublas_vector2.hpp>
#include<boost/numeric/bindings/traits/matrix_traits.hpp>
#include<boost/numeric/bindings/atlas/cblas.hpp>

#include"PGBasics.hh"
#include"NeuralNet.hh"
#include"NeuralNetBatch.hh"
#include"NeuralNetAtlas.hh"

namespace atlas = boost::numeric::bindings::atlas;

namespace libpg {

/**
 * Feedforward stage of the neural network. If you don't know how
 * neural nets work, don't even bother trying to understand this bit
 * of the code.
 * @param  input observation
 * @param  output approximation destination matrix
 */
void NeuralNetAtlas::doApprox(Observation& obs, Vector& output) {

    assert(obs.getFeatures().size1() == (unsigned int)dims(0));
    // For efficiency (both space and time) we treat the single layer
    // case and multi-layer cases differently, using the obs as the
    // first activation and output as the final output activation.

    //Linear network
    if (layers < 2) {
	//axpy_prod(column(obs.getFeatures(), obs.getAgent()), layerParams[0], output, true);
	atlas::gemv(CblasTrans, 1.,
		    layerParams[0],
		    column(obs.getFeatures(), obs.getAgent()),
		    0.,
		    output
		    );
	if ((int)squash[1]) squashVec(output);
    } 
    else {
	int l=0;

	// First layer
	//axpy_prod(column(obs.getFeatures(), obs.getAgent()), layerParams[l],  activations[l], true);
	atlas::gemv(CblasTrans, 1.,
		    layerParams[l],
		    column(obs.getFeatures(), obs.getAgent()),
		    0.,
		    activations[l]
		    );
	// Squash output of first layer? A little confusing since we
	// must check squash[l+1] to see if we squash activations[l]
	if ((int)squash[l+1]) squashVec(activations[l]);
	l++;
	
	// Middle layers
	for (; l < layers-1; l++) {
	    //axpy_prod(activations[l-1], layerParams[l], activations[l], true);
	    atlas::gemv(CblasTrans, 1.,
			layerParams[l],
			activations[l-1],
			0.,
			activations[l]
			);
	    if ((int)squash[l+1]) squashVec(activations[l]);
	}
	
	// Output layer
	//axpy_prod(activations[l-1], layerParams[l], output, true);
	atlas::gemv(CblasTrans, 1.,
		    layerParams[l],
		    activations[l-1],
		    0.,
		    output
		    );
	if ((int)squash[l+1]) squashVec(output);

    }

}


/**
 * Do error back propogation, but on the outputDeltas provided rather
 * than on any error MTL does not support matrix mult of 2 vectors to
 * a matrix directly (i.e., outer product), but does support this
 * through the rank_one_update command which seems quite efficient.
 * @param  input observation (again)
 * @param  output deltas, i.e., the gradient of the softmax function.
 * SIDEEFFECT: adds gradients to the trace variables in this class
 */
void NeuralNetAtlas::feedbackGrad(Observation& obs, Vector& outputDeltas) {
    
    // Single layer case
    if (layers < 2) {
	// rank one update = outer produc of two vectors
	// Compute grad for only layer
	//noalias(layerTraces[0]) += outer_prod(column(obs.getFeatures(), obs.getAgent()), outputDeltas);
	atlas::ger(column(obs.getFeatures(), obs.getAgent()),
		   outputDeltas, 
		   layerTraces[0]
		   );
    }
    else {
	int l = layers - 1;

	// Compute grad for last layer
	//noalias(layerTraces[l]) += outer_prod(activations[l-1], outputDeltas);
	atlas::ger(activations[l-1],
		   outputDeltas,
		   layerTraces[l]
		   );
	// Back prop deltas through last layer
	//axpy_prod(layerParams[l], outputDeltas, deltas[l-1], true);
	atlas::gemv(CblasNoTrans, 1.,
		    layerParams[l],
		    outputDeltas,
		    0.,
		    deltas[l-1]
		    );
	// Feedback through squashing
	if ((int)squash[l]) dSquashVec(activations[l-1], deltas[l-1]);
	l--;

	// repeat for all hidden layers
	for (; l > 0; l--) {
	    // Compute grad for current layer
	    //noalias(layerTraces[l]) += outer_prod(activations[l-1], deltas[l]);
	    atlas::ger(activations[l-1],
		       deltas[l-1],
		       layerTraces[l]
		       );
	    // Back prop deltas through current layer
	    //axpy_prod(layerParams[l], deltas[l], deltas[l-1], true);
	    atlas::gemv(CblasNoTrans, 1.,
			layerParams[l],
			deltas[l],
			0.,
			deltas[l-1]
			);
	    // Feedback through squashing
	    if ((int)squash[l]) dSquashVec(activations[l-1], deltas[l-1]);
	}

	// Compute grad for input layer
	//noalias(layerTraces[l]) += outer_prod(column(obs.getFeatures(), obs.getAgent()), deltas[l]);
	atlas::ger(column(obs.getFeatures(),obs.getAgent()),
		   deltas[l],
		   layerTraces[l]
		   );
	
    }

}

}
#endif

/**
 * Get rid of linking warning about no symbols when we don't have atlas.
 */
void stupidSymbolToShutUpGCCLinkerError(){}
