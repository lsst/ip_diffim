// -*- lsst-c++ -*-
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp> 
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/timer.hpp> 
#include <vw/Math/Matrix.h> 
#include "lsst/afw/Exception.h"

using namespace boost::numeric;
using namespace std;

int main(int argc, char** argv)
{
    int N = atoi(argv[1]);
    double e1, e2;

    /* just test matrices */
    ublas::matrix<double> b1(N,N), b2(N,N), b3(N,N);
    vw::Matrix<double>    v1(N,N), v2(N,N), v3(N,N);

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            e1 = (N * i + j) * rand();
            e2 = (N * i + j) * 2 * rand();

            b1(i, j) = e1;
            v1(i, j) = e1;

            b2(i, j) = e2;
            v2(i, j) = e2;
        }
    }
    

    boost::timer t;
    double time;

    // Sum
    t.restart();
    b3 = b1 + b2;
    time = t.elapsed();
    cout << "Boost sum " << time << "s" << endl; 

    t.restart();
    v3 = v1 + v2;
    time = t.elapsed();
    cout << "VW sum " << time << "s" << endl << endl; 
    

    // Multiplication
    t.restart();
    b3 = ublas::prod< ublas::matrix<double> >(b1,b2);
    time = t.elapsed();
    cout << "Boost multiplication " << time << "s" << endl; 

    t.restart();
    v3 = v1 * v2;
    time = t.elapsed();
    cout << "VW multiplication " << time << "s" << endl << endl; 


    // Transpose
    t.restart();
    b3 = ublas::trans< ublas::matrix<double> >(b1);
    time = t.elapsed();
    cout << "Boost transpose1 " << time << "s" << endl; 

    t.restart();
    b3 = ublas::trans< ublas::matrix<double> >(b2);
    time = t.elapsed();
    cout << "Boost transpose2 " << time << "s" << endl << endl; 

    t.restart();
    v3 = transpose(v1);
    time = t.elapsed();
    cout << "VW transpose1 " << time << "s" << endl; 

    t.restart();
    v3 = transpose(v2);
    time = t.elapsed();
    cout << "VW transpose2 " << time << "s" << endl << endl; 


    return 0;
}
