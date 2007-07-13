#include <vw/Math/Matrix.h> 
#include <PCA.h>
#include <iostream.h>

using namespace std;

int main( int argc, char** argv )
{
    
    typedef double ValueT;
    const int rowsN = 4;
    const int colsN = 5;
    vw::math::Matrix<ValueT> M(rowsN, colsN);
    // These matrices start at upper left
    M(0,0) = 1;
    M(0,4) = 2;
    M(1,2) = 3;
    M(3,1) = 4;
    cout << M << endl;

    vw::math::Matrix<ValueT> eVec(rowsN, colsN);
    vw::math::Vector<ValueT> eVal(colsN);

    lsst::imageproc::computePCAviaSVD(M, eVal, eVec);

    cout << "EVal " << eVal << endl;
    cout << "EVec " << vw::math::transpose(eVec) << endl;
    cout << "EVec " << vw::math::transpose(eVec)[0] << endl;
    exit(1);











    // screw it for now.  just work on it here.
//    vw::math::Matrix<double> u, vt;
//    vw::math::Vector<double> s;
//    vw::math::complete_svd( M, u, s, vt );
//
//    for (int i = 0; i < s.size(); i++) {
//        s[i] = vw::math::sqrt(s[i]);
//        cout << i << " " << s[i] << endl;
//    }
//
//    for (int col = 0; col < s.size(); col++) {
//        cout << "Eigenvector " << col << " : ";
//        for (int row = 0; row < s.size(); row++) {
//            cout << u(row, col) << " ";
//        }
//        cout << endl;
//    }
// 
//    boost::timer t;
//    double time;
//
//    // lets do some timing tests
//    int const nPixels = 30 * 30;
//    int const nSources = 30*30;
//    vw::math::Matrix<ValueT, nSources, nPixels> M2;
//    for (int i = 0; i < nSources; i++) {
//        for (int j = 0; j < nPixels; j++) {
//            M2(i, j)  = nPixels * i + j + 10 * rand();
//        }
//    }
//            
//
//    vw::math::Matrix<double> u2, vt2;
//    vw::math::Vector<double> s2;
//
//    t.restart(); 
//    //vw::math::complete_svd( M2, u2, s2, vt2 );
//    time = t.elapsed(); 
//
//    cout << "SVD took " << (boost::format("%.4f") % time) << "s" << endl; 

}
