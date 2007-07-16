#include <vw/Math/Matrix.h> 
#include <vw/Math/LinearAlgebra.h> 
#include <PCA.h>
#include <iostream>

using namespace std;

int main( int argc, char** argv )
{
    
    typedef double ValueT;
    typedef complex<ValueT> ValueTi;
    
    // Lets do a basic eigenvalue problem
    //    (0 1 0)
    // a =(1 0 0)
    //    (0 0 0)
    //
    // eigenvalues are -1, 0, 1
    //
    // eigenvectors are (1, -1, 0)
    //                  (0,  0, 1)
    //                  (1,  1, 0)
    const int rowsN = 3;
    const int colsN = 3;
    vw::math::Matrix<ValueT> M(rowsN, colsN);
    M(0, 0) = 0;
    M(0, 1) = 1;
    M(0, 2) = 0;
    
    M(1, 0) = 1;
    M(1, 1) = 0;
    M(1, 2) = 0;
    
    M(2, 0) = 0;
    M(2, 1) = 0;
    M(2, 2) = 0;
    cout << M << endl;
    
    cout << endl << "TESTING EIGENVECTOR CALCULATIONS" << endl;
    
    vw::math::Matrix<ValueTi> eVeci(rowsN, colsN);
    vw::math::Vector<ValueTi> eVali(colsN);
    vw::math::eigen(M, eVali, eVeci);
    
    // First test for the configuration of eigenvectors in eVeci
    for (int i = 0; i < eVeci.cols(); i++) {
        cout << "Eval " << i << " : " << eVali[i].real() << endl;
        cout << "Evec row " << i << " : " << vw::math::select_row(eVeci, i) << endl;
        cout << "Evec col " << i << " : " << vw::math::select_col(eVeci, i) << endl;
        cout << endl;
    }
    // Answers :
    cout << "ANSWER : Eval value -1 has Evec (+1, -1,  0)" << endl;
    cout << "ANSWER : Eval value  0 has Evec ( 0,  0, +1)" << endl;
    cout << "ANSWER : Eval value +1 has Evec (+1, +1,  0)" << endl;
    cout << "This means that the Eigenvectors are in the columns of Evec" << endl;
    cout << endl;
    
    // Compare SVD Results to Eigen
    cout << "Compare SVD Results to Eigen Results of M Mt" << endl;
    cout << " Tweaking M a bit" << endl;
    M(2, 1) = 7;
    M(1, 1) = 5;
    
    vw::math::Matrix<ValueT> MMt(rowsN, colsN);
    MMt = M * vw::math::transpose(M);
    
    vw::math::Matrix<ValueT> eVec(rowsN, colsN);
    vw::math::Vector<ValueT> eVal(colsN);
    vw::math::Vector<ValueT> colMean(colsN);
    lsst::imageproc::computePCA(MMt, colMean, eVal, eVec, false);
    
    vw::math::eigen(MMt, eVali, eVeci);
    
    // First test for the configuration of eigenvectors in eVeci
    for (int i = 0; i < eVeci.cols(); i++) {
        cout << "Eval " << i << " : " << endl;
        cout << " SVD     : " << eVal[i] << endl;
        cout << " Eig     : " << eVali[i].real() << endl;
        
        cout << "Evec " << i << " : " << endl;
        
        cout << " SVD row : " << vw::math::select_row(eVec, i) << endl;
        cout << " Eig row : " << vw::math::select_row(eVeci, i) << endl;
        cout << " SVD col : " << vw::math::select_col(eVec, i) << endl;
        cout << " Eig col : " << vw::math::select_col(eVeci, i) << endl;
        
        cout << endl;
    }
    
    cout << "CONCLUSION : SVD eigenvectors are in the columns of eVec, and are sorted by Eigenvalue and normalized" << endl << endl;
    
}
