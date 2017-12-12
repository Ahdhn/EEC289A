#include"PGBasics.hh"
#include"Teacher.hh"
#include"Bias.hh"

#define OutputWeight (4.0)

using namespace std;
namespace libpg {

    /**
     * Use a teacher to add an amount of bias to an approximators output.
     */
    Bias::Bias(Approximator* approx , Teacher* teacher) : TransformApproximator(approx) {

	/*assert(teacher->getOutputDim() == approx->getOutputDim());
	 *
	 * Note: at this stage, the teacher may not be completely initialised.
	 * In particular, its output dimension may not be set properly.
	 */
	this->teacher = teacher;
	teacherOutput.resize(approx->getOutputDim());

    }


    /**
     * Where the real coding is done to implement softmax
     */
    void Bias::doApprox(Observation& obs, Vector& output) {

	approx->doApprox(obs, output);
	teacher->doApprox(obs, teacherOutput);
	output += OutputWeight * teacherOutput;
    }
    
    
}
