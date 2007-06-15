// -*- lsst-c++ -*-
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>

using namespace boost::numeric;
using namespace std;

int main()
{
    typedef ublas::vector<double> Vector;
    typedef ublas::matrix<double> Matrix;

    /* just a stupid test matrix */
    Matrix m(3,3);
    Matrix inverse(3,3);

    for (unsigned i = 0; i < 3; ++ i)
        for (unsigned j = 0; j < 3; ++ j)
            m (i, j) = 3 * i + j + 2;
    
    //Matrix I3 = ublas::identity_matrix<double>(3, 3);
    //m = m + 3 * I3; // now m is not singular
    
    // create a working copy of the input
    Matrix<double> A(m);
    
    // create a permutation matrix for the LU-factorization
    ublas::permutation_matrix<std::size_t> pm(A.size1());

    // perform LU-factorization
    int res = lu_factorize(A, pm);
    if( res != 0 ) return false;
    
    // create identity matrix of "inverse"
    inverse.assign(ublas::identity_matrix<double>(A.size1()));
    
    // backsubstitute to get the inverse
    lu_substitute(A, pm, inverse);
    
    // can i do this?
    Matrix p(m * pm);
    cout << " input = " << m << endl;
    cout << " inverse = " << pm << endl;
    cout << " product = " << p << endl;
    
    return 0;
} 
