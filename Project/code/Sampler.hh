#ifndef Sampler_hh
#define Sampler_hh
// $Id: Sampler.hh 127 2007-09-10 17:20:04Z daa $

namespace libpg {
class Sampler {

public:
    
    static int discrete(Vector& pdf, double cdf);

}; 
}
#endif
