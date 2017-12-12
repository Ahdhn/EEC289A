#include"PGBasics.hh"
#include"Sampler.hh"

using namespace std;

namespace libpg {

int Sampler::discrete(Vector& pdf, double cdf) {

    double thresh = random()/(double)RAND_MAX*cdf;
    double partialCDF = 0.0;

    Vector::iterator i;
    for (i = pdf.begin(); i != pdf.end(); i++) {
	partialCDF += *i;
	if (partialCDF >= thresh) break; 
    }

    if (i == pdf.end()) {
	// Should never reach the end. Although sometimes happens if 
	// approximator output goes NaN
	cerr<<"!! Sampler broke: pdf="<<pdf<<" cdf="<<cdf<<" thresh="<<thresh<<" partialCDF="<<partialCDF<<endl; 
	free((void*)1);
	throw std::runtime_error("Sampler broke\n");
    }
    assert(i.index() < pdf.size());
    return i.index();

}
}
