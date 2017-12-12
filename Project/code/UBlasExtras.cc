// $Id: UBlasExtras.cc 127 2007-09-10 17:20:04Z daa $

#include"PGBasics.hh"
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/blas.hpp>

using namespace std;


#ifdef HAVE_BOOST_SANDBOX
/**
 * We can include some funky features like SVD, and perhaps Cholesky
 * least squares.
 */
#include <boost/numeric/bindings/traits/std_vector.hpp>
#include <boost/numeric/bindings/traits/type_traits.hpp>
#include <boost/numeric/bindings/traits/vector_traits.hpp>
#include <boost/numeric/bindings/traits/matrix_traits.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/traits/ublas_banded.hpp>
#include <boost/numeric/bindings/traits/ublas_hermitian.hpp>
#include <boost/numeric/bindings/traits/detail/symm_herm_traits.hpp>
#include <boost/numeric/bindings/traits/ublas_symmetric.hpp>
#include <boost/numeric/bindings/traits/ublas_vector.hpp>
#include <boost/numeric/bindings/traits/traits.hpp>
#include <boost/numeric/bindings/lapack/gesvd.hpp>
#include <boost/numeric/bindings/lapack/posv.hpp>
#include <boost/numeric/bindings/lapack/lapack.h>

#ifdef HAVECLAPACK
#include <boost/numeric/bindings/atlas/clapack.hpp>
#endif

#endif

namespace libpg {

/**
 * Create a random matrix with entrie
 * @param matrix to randomize
 * @param maximum absolute value. Negative entries are also made.
 */ 
void UBlasExtras::randomize(Matrix& m, double maxRand) {

    for (unsigned int i = 0; i < m.size1(); i++) {
	for (unsigned int j = 0; j < m.size2(); j++) {
	    m(i,j) = ((random()/(double)RAND_MAX) - 0.5)*2*maxRand;
	}
    }
    
}



/**
 * Matrix inversion routine.
 * Uses lu_factorize and lu_substitute in uBLAS to invert a matrix 
 * @param input matrix
 * @param destination for inverse.
 */
void UBlasExtras::invertMatrix(const Matrix& input, Matrix& inverse) {
    
    typedef permutation_matrix<std::size_t> pmatrix;
     // create a working copy of the input
     Matrix A(input);
     // create a permutation matrix for the LU-factorization
     pmatrix pm(A.size1());
     
     // perform LU-factorization
     lu_factorize(A,pm);
     
     // create identity matrix of "inverse"
     inverse.assign(identity_matrix<double>(A.size1()));
     
     // backsubstitute to get the inverse
     lu_substitute(A, pm, inverse);
 }



/**
 * Solve a linear system Ax=b for x. Assumes full rank.
 * @param A matrix
 * @param b vector
 * @param destination vector
 */
void UBlasExtras::solveLinear(Matrix& A, Vector& b, Vector& x) {

    typedef permutation_matrix<std::size_t> pmatrix;
    // create a working copy of the input
    Matrix input(A);
    // create a permutation matrix for the LU-factorization
    pmatrix pm(A.size1());
    
    // perform LU-factorization
    lu_factorize(input,pm);
    x = b;
    lu_substitute(input, pm, x);
}
    


/**
 *Use Sherman-Morrison fomular to update inverse of a Matrix
 *  Formula: (A+cd')^-1 = A^-1 - (A^-1cd'A^-1) / (1+d'A^-1c)
 * @param A^-1 from  previous step
 * @param vector c as in A_{t+1} = A_t + cd^T
 * @param vector d as in A_{t+1} = A_t + cd^T
 */
void UBlasExtras::updateInverseMatrix(Matrix& inverse, Vector& c, Vector& d) {

    Vector a(c.size());
    a.clear();
    Vector b(d.size());
    b.clear();
    axpy_prod(inverse,c,a,true);
    axpy_prod(d,inverse,b,true);
    inverse = inverse - outer_prod(a,b)/(1.0 + inner_prod(b,c));
    
}


/**
 * @param  matrix to convert to a vector 
 * @param  output will be the sum of this input plus vectorised matrix.
 * @param  requires, second arg dimension should be > size1()*size2()
 * of first input matrix
 * @param  starting point in the vector, this allows vectors to
 * aggregate many small matricies. Defaults to 0
 */
void UBlasExtras::addMatrixToVector(Matrix& m, Vector& v, int s) {

     assert(m.size1()*m.size2() <= v.size());
     
     for (size_t j=0; j<m.size2(); j++) {
	 project(v, range(j*m.size1() + s, (j+1)*m.size1() + s)) += column(m, j);
     }

 }


/**
 * split up scaled vector and add to matrix
 * INPUT scale factor
 * @param  output will be the sum of this input plus vectorised matrix.
 * @param  requires, second arg dimension should be > size1()*size2()
 * of first input matrix
 * @param  starting point in the vector, this allows vectors to
 * aggregate many small matricies. Defaults to 0
 */
void UBlasExtras::addScaledVectorToMatrix(double scale, Vector& v, Matrix& m, int s) {

     assert(m.size1()*m.size2() <= v.size());
     
     for (size_t j=0; j<m.size2(); j++) {
	 column(m, j) +=  scale*project(v, range(j*m.size1() + s, (j+1)*m.size1() + s));
     }

 }


/**
 * Compute the angle between two vectors
 * Not very efficient.
 * Output: radians
 */
double UBlasExtras::angle(Vector& v1, Vector& v2) {

    return acos(inner_prod(v1, v2)/(norm_2(v1)*norm_2(v2)));

}


/**
 *  Check in a vector for NaNs
 *  Pretty slow!
 */
bool UBlasExtras::containsNaN(Vector& v) {
    for (size_t i = 0; i < v.size(); i++) if (isnan(v[i])) return true;
    return false;
}


/**
 *  Check in a vector for NaNs
 *  Pretty slow!
 */
bool UBlasExtras::containsNaN(Matrix& m) {
    for (size_t i = 0; i < m.size1(); i++) {
        for (size_t j = 0; j < m.size2(); j++) {
	    if (isnan(m(i,j))) return true;
	}
    }
    return false;
}


/**
 * Perform least squares regression.
 * Should use cholesky solve here because A trans(A) is symmetric pos. 
 * @param Matrix A with one sample per row
 * @param Vector y with observed outcomes
 * @param destination vector
 */
void UBlasExtras::linearLeastSquares(Matrix& A, Vector& y, Vector& x) {

    Matrix sqrA(A.size2(), A.size2());
    Vector trAy(x.size()); 

    // Symmetry optimised rank k update.
    noalias(sqrA) = prod(trans(A), A);

    //cout<<"sqrA:"<<sqrA<<endl;

    noalias(trAy) = prod(trans(A), y);
    //cout<<"trAy:"<<trAy<<endl;

#if (HAVE_BOOST_SANDBOX && HAVE_ATLAS) 
    solveLinearCholesky(sqrA, trAy, x);
#else
    solveLinear(sqrA, trAy, x);
#endif


}


void UBlasExtras::linearLeastSquares_test() {

    Matrix A(2,2);
    A.clear();
    A(0,0) = 1.0;
    A(1,1) = 1.0;

    Vector b(2);
    b[0] = 1.0;
    b[1] = 1.0;

    Vector x(2);
    
    linearLeastSquares(A, b, x);

    cout<<"linearLeastSquares_test() 1:"<<endl
	<<"\tA:"<<A<<endl
	<<"\tb:"<<b<<endl
	<<"\tx:"<<x<<endl;


    A.resize(4, 2);
    A(0,0) = 0; A(0, 1) = 1;
    A(1,0) = 2; A(1, 1) = 1;
    A(2,0) = 4; A(2, 1) = 1;
    A(3,0) = -1; A(3, 1) = 1;

    b.resize(4);
    b[0] = 3;
    b[1] = 3;
    b[2] = 4;
    b[3] = 2;

    x.resize(2);
    linearLeastSquares(A, b, x);

    cout<<"linearLeastSquares_test() 2:"<<endl
	<<"\tA:"<<A<<endl
	<<"\tb:"<<b<<endl
	<<"\tx:"<<x<<endl;


}


/**
 * Perform singular value demomposition using clapack and the ublas bindings.
 * @param vector to put singular values in.
 */
void UBlasExtras::svd(Matrix& A, Vector& s) {

#if (HAVE_BOOST_SANDBOX && HAVE_CLAPACK)


    int const m = boost::numeric::bindings::traits::matrix_size1(A);
    int const n = boost::numeric::bindings::traits::matrix_size2(A);
    int const min_mn = std::min(m, n);
    
    // temporary storage
    Matrix at(A);
    Matrix u(m, min_mn); 
    u.clear();
    Matrix vt(min_mn, n); 
    vt.clear();
    s.resize(min_mn);
    
    // We do have column major matricies.


    // call SVD
    int info = boost::numeric::bindings::lapack::gesvd(at, s, u, vt);
    if (info != 0) throw std::runtime_error("failed in gesvd");
#else
    throw std::runtime_error("SVD needs CLAPACK and BOOST_SANDBOX. Check makefiles.");
#endif

}



/**
 * Compute rank of matrix based on SVD.  Does not cope with noise
 * beyond simple threshold for errors at machine precision level.
 * @param Matrix to take SVD of. Rank will be less or equal to smallest diim
 * @returns rank >= 0, <= size of smallest dim of A.
 */
size_t UBlasExtras::rank(Matrix& A) {

    Vector sv;
    size_t r=0;
    size_t s=0;

    svd(A, sv);

    while (s < sv.size() && sv[s++] > 1e-15) r++;

    return r;
	
}


/**
 * Compute pseudo inverse from SVD.
 * If A = USB', where S is diagonal, then pseudo inverse is
 * Ainv = BS^{-1}A'
 */ 
Matrix& UBlasExtras::pseudoInverse(Matrix& A, Matrix& Ainv) {

#if (HAVE_BOOST_SANDBOX && HAVE_CLAPACK)


    int const m = boost::numeric::bindings::traits::matrix_size1(A);
    int const n = boost::numeric::bindings::traits::matrix_size2(A);
    int const min_mn = std::min(m, n);
    
    // temporary storage
    Matrix at(A);
    Matrix u(m, min_mn); 
    u.clear();
    Matrix vt(min_mn, n); 
    vt.clear();
    Vector s(min_mn);
    s.clear();

    // We do have column major matricies.


    // call SVD
    int info = boost::numeric::bindings::lapack::gesvd(at, s, u, vt);
    if (info != 0) throw std::runtime_error("failed in gesvd");

    at.resize(min_mn, min_mn);
    at.clear();
    for (size_t sv=0; sv < s.size(); sv++) {
	if (s[sv] > 0.0) at(sv, sv) = 1.0/s[sv];
	else at(sv, sv)=0;
    }

    assert(Ainv.size1() == A.size2());
    assert(Ainv.size2() == A.size1());
    Matrix tmp(min_mn, m);
    noalias(tmp) = prod(at, trans(u));
    Ainv = prod(trans(vt), tmp);


#else
    throw std::runtime_error("SVD needs CLAPACK and BOOST_SANDBOX. Check makefiles.");
#endif

    return Ainv;

}



/**
 * Perform singular value demomposition using clapack and the ublas bindings.
 * @param vector to put singular values in.
 */
void UBlasExtras::solveLinearCholesky(Matrix& A, Vector& b, Vector& x) {

#if (HAVE_BOOST_SANDBOX && HAVE_CLAPACK) 
    // temporary storage

    symmetric_adaptor<Matrix, lower> symA(A);
    Matrix bt(b.size(), 1);
    column(bt, 0).assign(b);
    // We do have column major matricies.

    // call Cholesky solve. Not sure why, but on OSX clapack (with
    // ATLAS BLAS) seems to beat ATLAS.
    int info = boost::numeric::bindings::lapack::posv(symA, bt);
    if (info != 0) throw std::runtime_error("failed in posv, Cholesky solve");

    noalias(x) = column(bt, 0);

#else
    throw std::runtime_error("Cholesky solve needs CLAPACK. Check makefiles for -DHAVE_CLAPACK.");
#endif

}
}
