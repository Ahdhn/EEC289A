// $Id$

#include "PGBasics.hh"
#include "TransformApproximator.hh"

namespace libpg {


    /**
     * Create a wrapper around an existing approximator only.
     * @param  approx base approximator
     */
    TransformApproximator::TransformApproximator(Approximator* approx) {
	this->approx = approx;
    }

        
    /** 
     * Every single method here just calls its base equivalent.
     */
    int TransformApproximator::getNumParams() { return approx->getNumParams(); }
    void TransformApproximator::doApprox(Observation& obs, Vector& output) { approx->doApprox(obs, output); }
    void TransformApproximator::feedbackGrad(Observation& obs, Vector& deltas) { approx->feedbackGrad(obs, deltas); }
    void TransformApproximator::discountTrace() { approx->discountTrace(); };
    void TransformApproximator::setDiscount(double d) { approx->setDiscount(d); }
    void TransformApproximator::instantStep(double d) { approx->instantStep(d); }
    void TransformApproximator::setStepSize(double d) { approx->setStepSize(d); }
    void TransformApproximator::resetTrace() { approx->resetTrace(); }
    void TransformApproximator::resetParams() { approx->resetParams(); }
    void TransformApproximator::randomizeParams(double d) {approx->randomizeParams(d); }
    double TransformApproximator::getMaxParam() { return approx->getMaxParam(); }
    void TransformApproximator::write(std::ostream& os) { approx->write(os); }
    void TransformApproximator::read(std::istream& is) { approx->read(is); }
    int TransformApproximator::getInputDim() { return approx->getInputDim(); }
    int TransformApproximator::getOutputDim() { return approx->getOutputDim(); }
    void TransformApproximator::batchStep() { approx->batchStep(); }
    void TransformApproximator::accumulateGrad(double d, Observation& obs) { approx->accumulateGrad(d, obs); }
    void TransformApproximator::resetGrad() { approx->resetGrad(); }
    void TransformApproximator::computeDirection(int i) { approx->computeDirection(i); }
    void TransformApproximator::reduce(Vector& v, Approximator::StatsEnum s) { approx->reduce(v, s); }
    void TransformApproximator::scatter(Vector& v, Approximator::StatsEnum s) { approx->scatter(v, s); }

}
