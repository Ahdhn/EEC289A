#ifndef UBlasExtras_hh
#define UBlasExtras_hh
// $Id: UBlasExtras.hh 127 2007-09-10 17:20:04Z daa $
namespace libpg {

/**
 * Static methods to augment the general UBLAS library methods.
 */

class UBlasExtras {

public:

    static bool containsNaN(Matrix& m);
    static bool containsNaN(Vector& m);
    static void randomize(Matrix& m, double maxRand);
    static void invertMatrix(const Matrix& input, Matrix& inverse);
    static void solveLinear(Matrix& A, Vector& b, Vector& x);
    static void solveLinearCholesky(Matrix& A, Vector& b, Vector& x);
    static void updateInverseMatrix(Matrix& inverse, Vector& c, Vector& d); 
    static void addMatrixToVector(Matrix& m, Vector& v, int s=0);
    static void addScaledVectorToMatrix(double scale, Vector& v, Matrix& m, int s=0);
    static double angle(Vector& v1, Vector& v2);
    static void linearLeastSquares(Matrix& A, Vector& y, Vector& x);
    static void svd(Matrix& A, Vector& s);
    static size_t rank(Matrix& A);
    static Matrix& pseudoInverse(Matrix& A, Matrix& Ainv);

    static void linearLeastSquares_test();
};
}
#endif
