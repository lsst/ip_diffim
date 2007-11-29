// -*- lsst-c++ -*-
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp> 
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/timer.hpp> 
#include <vw/Math/Matrix.h> 
#include <vw/Math/LinearAlgebra.h> 
#include <lsst/mwi/exceptions/Exception.h>

//using namespace boost::numeric;  // Watch namespace issues; make sure "invert" is VW and not boost
using namespace std;

int main(int argc, char** argv)
{
    int N = atoi(argv[1]);

    /* just a stupid test matrix */
    boost::numeric::ublas::matrix<double> m(N,N), inv(N,N), p(N,N);
    vw::Matrix<double> mvw(N,N);
    
    
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            m (i, j)  = N * i + j + 2;
            mvw(i, j) = N * i + j + 2;
        }
    }
    
    boost::numeric::ublas::matrix<double> I = boost::numeric::ublas::identity_matrix<double>(N, N);
    m = m + N * I; // now m is not singular

    vw::Matrix<double> Ivw(N,N);
    Ivw.set_identity();
    mvw = mvw + N * Ivw; 
    
    boost::timer t;
    double time;

    // BELOW IS THE GUTS OF THE CODE 
    //
    //
    t.restart(); 

    // create a working copy of the input
    boost::numeric::ublas::matrix<double> A(m);
    
    // create a permutation matrix for the LU-factorization; pivots
    boost::numeric::ublas::permutation_matrix<std::size_t> pm(A.size1());

    // perform LU-factorization; for LSST throw exception instead of return
    /*
    try {
        int res = lu_factorize(A, pm);
        if( res != 0 ) throw lsst::mwi::exceptions::Runtime("LU factorization failed; singular matrix.  Get a better exception than this one!");
    } catch (lsst::mwi::exceptions::Exceptionstack &e) {
        cerr << "Caught exception: " << e.what() << "\n";
        return 1;
    }
    */

    // create identity matrix of "inverse"
    inv.assign(boost::numeric::ublas::identity_matrix<double>(A.size1()));
    
    // backsubstitute to get the inverse
    lu_substitute(A, pm, inv);

    time = t.elapsed(); 
    cout << "Boost took " << time << "s" << endl; 
    //
    //
    // ABOVE IS THE GUTS OF THE CODE 

    // verify you get out the identity matrix
    p = boost::numeric::ublas::prod< boost::numeric::ublas::matrix<double> >(m,inv); 

    if (N < 5) {
        cout << " input = " << m << endl;
        cout << " inverse = " << inv << endl;
        cout << " product = " << p << endl;
    }

    // TRY THE SAME THING WITH VW
    // 
    // 
    t.restart(); 
    vw::Matrix<double> Avw = vw::math::inverse(mvw);
    time = t.elapsed(); 
    cout << "VW took " << time << "s" << endl;     
    if (N < 5) {
        cout << " input = " << mvw << endl;
        cout << " inverse = " << Avw << endl;
    }

    // TRY THE SAME THING WITH VW
    // 
    // 
    t.restart(); 
    vw::Matrix<double> Avw2 = vw::math::pseudoinverse(mvw);
    time = t.elapsed(); 
    cout << "VW took " << time << "s" << endl;     
    if (N < 5) {
        cout << " input = " << mvw << endl;
        cout << " inverse = " << Avw2 << endl;
    }

    return 0;
} 
