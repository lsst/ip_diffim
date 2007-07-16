#include <vw/Math/Matrix.h> 
#include <vw/Math/LinearAlgebra.h> 
#include <PCA.h>
#include <iostream>

using namespace std;

int main( int argc, char** argv )
{
    
    typedef double ValueT;
    // Lets try decomposing another matrix using Basis functions
    const int rowsN2 = 4; // Number of variables = m
    const int colsN2 = 5; // Number of realizations = n
    vw::math::Matrix<ValueT> M2(rowsN2, colsN2); // m x n
    vw::math::Matrix<ValueT> eVec2(rowsN2, colsN2); // m x n
    vw::math::Vector<ValueT> eVal2(colsN2); // n 
    vw::math::Vector<ValueT> rowMean2(rowsN2); // m
    M2(0,0) = 1;
    M2(0,4) = 2;
    M2(1,2) = 3;
    M2(3,1) = 4;
    cout << M2 << endl;

    lsst::imageproc::computePCA(M2, rowMean2, eVal2, eVec2);
    
    cout << endl << "NEW PCA" << endl;
    for (int i = 0; i < eVec2.cols(); i++) {
        cout << "Eval " << i << " : " << eVal2[i] << endl;
        cout << "Evec col " << i << " : " << vw::math::select_col(eVec2, i) << endl;
        cout << endl;
    }
    
    int nCoeff = 4;
    vw::math::Matrix<ValueT> coeff2(colsN2, rowsN2);

    lsst::imageproc::decomposeMatrixUsingBasis(M2, eVec2, coeff2);
    cout << "C   : " << coeff2 << endl;
    cout << "In  : " << M2 << endl;
    vw::math::Matrix<ValueT> M3(rowsN2, colsN2); // m x n
    lsst::imageproc::approximateMatrixUsingBasis(eVec2, coeff2, nCoeff, M3);
    cout << "Out : " << M3 << endl;

    cout << endl;
    
    vw::math::Matrix<ValueT> coeff3(colsN2, nCoeff);
    lsst::imageproc::decomposeMatrixUsingBasis(M2, eVec2, nCoeff, coeff3);
    cout << "C   : " << coeff3 << endl;
    lsst::imageproc::approximateMatrixUsingBasis(eVec2, coeff3, nCoeff, M3);
    cout << "Out : " << M3 << endl;

}











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


