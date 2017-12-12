#ifndef PGBasics_hh
#define PGBasics_hh


/**
 * \mainpage The PG Library 
 * 
 * \verbinclude README.txt
 */


/**
 * Set up the Matrix Template Library with a basic matrix and vector
 */
#include <iostream>
#include <fstream>
#include <assert.h>
#include <vector>
#include <list>
#include <set>
#include <map>
#include <string>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/operation.hpp>

// Set up the default name spaces

using namespace boost::numeric::ublas;


namespace libpg {

/**
 * Normal matrix
 */
typedef matrix<double, column_major> Matrix;

/**
 * Symmetric matrix, mostly for least squares stuff.
 */
typedef symmetric_matrix<double, lower, column_major> SymmetricMatrix;

/**
 * Mapped Matrix for sparse data.
 */
typedef mapped_matrix<double, column_major> SparseMatrix;
 
/**
 * Normal vector
 */
typedef boost::numeric::ublas::vector<double> Vector;

/**
 * Boolean vector
 */
typedef boost::numeric::ublas::vector<bool> bVector;
    //typedef std::vector<bool> bVector; // <- DON'T USE THIS ! THIS IS SLOW...

/**
 * Diagonal matrix
 */
typedef banded_matrix<double> Diagonal; 

// Some handy macros
#ifndef NYI
#define NYI { std::cerr<<__FUNCTION__<<" is NYI in "<<__FILE__<<":"<<__LINE__<<std::endl; abort(); }
#endif

#define CWIDTH 12 // Column width for output
#define PRECISION 5 // Precision for output

#define PG_MACHINE_EPS 1e-10

}

// Include basic header files we will need for any PG app
#include "UBlasExtras.hh"
#include "Observation.hh"
 
#include "Approximator.hh"
#include "Policy.hh"
#include "Controller.hh"
#include "TransformController.hh"
#include "Simulator.hh"
#include "RLAlg.hh"

#endif

