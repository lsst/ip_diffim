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
    const int colsN2 = 6; // Number of realizations = n
    vw::math::Matrix<ValueT> M2(rowsN2, colsN2); // m x n
    vw::math::Matrix<ValueT> eVec2(rowsN2, colsN2); // m x n
    vw::math::Vector<ValueT> eVal2(colsN2); // n 
    vw::math::Vector<ValueT> rowMean2(rowsN2); // m
    M2(0,0) = 1;
    M2(0,4) = 2;
    M2(1,2) = 3;
    M2(3,1) = 4;
    M2(2,5) = 6;

    // Test case where colsN2 > rowsN2
    // 
    // 
    cout << endl << "COLS > ROWS" << endl;
    cout << M2 << endl;

    lsst::imageproc::computePCA(M2, rowMean2, eVal2, eVec2);
    
    for (unsigned int i = 0; i < eVec2.cols(); i++) {
        cout << "Eval " << i << " : " << eVal2[i] << endl;
        cout << "Evec col " << i << " : " << vw::math::select_col(eVec2, i) << endl;
        cout << endl;
    }

    vw::math::Matrix<ValueT> coeff2(colsN2, colsN2); // n x n
    vw::math::Matrix<ValueT> M2a(rowsN2, colsN2); // m x n

    cout << "Input data to reconstruct : " << endl;
    cout << M2 << endl;

    int nCoeff;
    if (colsN2 < rowsN2) {
        nCoeff = colsN2;
    }
    else {
        nCoeff = rowsN2;
    }

    // test with Ncoeff
    for (int i = nCoeff; i > 0; i--) {
        coeff2.set_size(colsN2, i);

        cout << " N coeff = " << i << endl;
        lsst::imageproc::decomposeMatrixUsingBasis(M2, eVec2, i, coeff2);
        cout << " Coeff : " << endl;
        cout << "  " << coeff2 << endl;
        lsst::imageproc::approximateMatrixUsingBasis(eVec2, coeff2, i, M2a);
        cout << " Out   : " << endl;
        cout << "  " << M2a << endl;
        cout << " Diff  : " << endl;
        cout << "  " << M2-M2a << endl << endl;
    }


    // Test case where rowsN2 > colsN2
    // 
    // 
    int rowsN3 = colsN2;
    int colsN3 = rowsN2;
    vw::math::Matrix<ValueT> M3(rowsN3, colsN3); // m x n
    vw::math::Matrix<ValueT> eVec3(rowsN3, colsN3); // m x n
    vw::math::Vector<ValueT> eVal3(colsN3); // n 
    vw::math::Vector<ValueT> rowMean3(rowsN3); // m
    M3(0,0) = 1;
    M3(4,0) = 2;
    M3(2,1) = 3;
    M3(1,3) = 4;
    M3(5,2) = 6;
    cout << endl << "ROWS > COLS" << endl;
    cout << M3 << endl;

    lsst::imageproc::computePCA(M3, rowMean3, eVal3, eVec3);
    
    for (unsigned int i = 0; i < eVec3.cols(); i++) {
        cout << "Eval " << i << " : " << eVal3[i] << endl;
        cout << "Evec col " << i << " : " << vw::math::select_col(eVec3, i) << endl;
        cout << endl;
    }

    vw::math::Matrix<ValueT> coeff3(colsN3, colsN3); // n x n
    vw::math::Matrix<ValueT> M3a(rowsN3, colsN3); // m x n

    cout << "Input data to reconstruct : " << endl;
    cout << M3 << endl;

    if (colsN3 < rowsN3) {
        nCoeff = colsN3;
    }
    else {
        nCoeff = rowsN3;
    }

    // test with Ncoeff
    for (int i = nCoeff; i > 0; i--) {
        coeff3.set_size(colsN3, i);

        cout << " N coeff = " << i << endl;
        lsst::imageproc::decomposeMatrixUsingBasis(M3, eVec3, i, coeff3);
        cout << " Coeff : " << endl;
        cout << "  " << coeff3 << endl;
        lsst::imageproc::approximateMatrixUsingBasis(eVec3, coeff3, i, M3a);
        cout << " Out   : " << endl;
        cout << "  " << M3a << endl;
        cout << " Diff  : " << endl;
        cout << "  " << M3-M3a << endl << endl;
    }
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


