// -*- lsst-c++ -*-
/**
 * \file
 *
 * Implementation of PCA Methods
 *
 * Implementation of PCA Methods
 *
 * \author Andrew Becker
 *
 * \ingroup imageproc
 */

#include <vw/Math/Matrix.h> 
#include <vw/Math/Vector.h> 
#include <vw/Math/LinearAlgebra.h> 
#include <vw/Math/Functions.h> 
#include <iostream>
using namespace std;

template <typename aMatrixT, typename aVectorT>
void lsst::imageproc::computePCAviaSVD(
    aMatrixT &M,
    aVectorT &eVal,
    aMatrixT &eVec
    ) {

    // Use LAPACK SVD
    // http://www.netlib.org/lapack/lug/node32.html
    // Suggests to me that we use DGESVD or DGESDD (preferred)
    // Good, gesdd is used in LinearAlgebra.h

    // All computations here are in double
    // Cast to aMatrixT and aVectorT after computation
    // This might be unncessarily inefficient
    vw::math::Matrix<double> u, vt;
    vw::math::Vector<double> s;
    vw::math::complete_svd(M, u, s, vt);

    /* Note on the math :

    In the SVD, we find matrices U S and V such that 
       M = U S V*
    where V* is the conjugate transpose of V.  For a real-valued matrix this is just Vt.  
    The diagonal entries of S are the singular values of M.

    In this decomposition
       M* M = V S* U* U S V*
            = V (S* S) V*

       M M* = U S V* V S* U*
            = U (S S*) U*

    ----------

    In an eigenvalue decomposition, we assume that
       M = U L U*
    where U is a Unitary matrix, and L is the diagonal matrix with the eigenvalues of M.

    ----------

    By comparing these results, we can see that

     o The SVD columns of U represent the eigenvectors of M M*
     o The SVD columns of V represent the eigenvectors of M* M
     o The SVD diagonal entries of S are the square root of the eigenvalues of both M M* and M* M

    ----------

    Finally, in a PCA, we want to find the eigenvectors of the covariance matrix M M*
    Therefore the SVD yields principal components with the eigenvectors in U

    */

    // Have s represent the eigenvalues; they are already sorted by LAPACK
    for (int i = 0; i < s.size(); i++) {
        eVal[i] = vw::sqrt(s[i]);
    }
    // We could use VectorProxys to do this

    // Eigenvectors are in the columns of eVec
    for (int col = 0; col < u.cols(); col++) {
        for (int row = 0; row < u.rows(); row++) {
            eVec(row,col) = u(row, col);
        }
    }
    // We could use MatrixProxys to do this 
}

template <typename aMatrixT>
void lsst::imageproc::decomposeMatrixUsingBasis(
    aMatrixT &M,
    aMatrixT &eVec,
    aMatrixT &coeff
    ) {
    
    // We assume that the *rows* of M hold the input signal
    // And the *cols* of eVec hold the basis functions

    //coeff.set_size(M.rows(), eVec.cols());
    
    for (int row = 0; row < M.rows(); row++) {
        vw::math::Vector<double> mRow = vw::math::select_row(M, row);
        cout << "Row " << row << " " << mRow << endl;
        for (int col = 0; col < eVec.cols(); col++ ) {
            vw::math::Vector<double> eCol = vw::math::select_col(eVec, col);
            
            cout << "Evec " << col << " " << eCol << endl;

            coeff(row, col) = vw::math::dot_prod(mRow, eCol);

            cout << "  Prod " << vw::math::dot_prod(mRow, eCol) << endl;

        }
    }
}
