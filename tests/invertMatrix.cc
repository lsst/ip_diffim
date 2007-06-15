// -*- lsst-c++ -*-
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp> 
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include "lsst/fw/Exception.h"

using namespace boost::numeric;
using namespace std;

int main()
{
    typedef ublas::vector<double> Vector;
    typedef ublas::matrix<double> Matrix;

    /* just a stupid test matrix */
    Matrix m(3,3), inverse(3,3), p(3,3);

    for (unsigned i = 0; i < 3; ++ i)
        for (unsigned j = 0; j < 3; ++ j)
            m (i, j) = 3 * i + j + 2;
    
    Matrix I3 = ublas::identity_matrix<double>(3, 3);
    m = m + 3 * I3; // now m is not singular
    
    // BELOW IS THE GUTS OF THE CODE 
    //
    //
    // create a working copy of the input
    Matrix A(m);
    
    // create a permutation matrix for the LU-factorization
    ublas::permutation_matrix<std::size_t> pm(A.size1());

    // perform LU-factorization; for LSST throw exception instead of return
    try {
        int res = lu_factorize(A, pm);
        if( res != 0 ) throw lsst::fw::Memory("LU factorization failed; singular matrix.  Get a better exception than this one!");
    } catch (lsst::fw::Exception &e) {
        cerr << "Caught exception: " << e.what() << "\n";
        return 1;
    }

    // create identity matrix of "inverse"
    inverse.assign(ublas::identity_matrix<double>(A.size1()));
    
    // backsubstitute to get the inverse
    lu_substitute(A, pm, inverse);
    
    //
    //
    // ABOVE IS THE GUTS OF THE CODE 

    // verify you get out the identity matrix
    p = ublas::prod< ublas::matrix<double> >(m,inverse); 

    cout << " input = " << m << endl;
    cout << " inverse = " << inverse << endl;
    cout << " product = " << p << endl;
    
    return 0;
} 
