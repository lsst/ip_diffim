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

using namespace std;

template <typename aMatrixT, typename aVectorT>
void lsst::imageproc::computePCAviaSVD(
    aMatrixT &A,
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
    vw::math::complete_svd(A, u, s, vt);

    /* Note on the math :

    In the SVD, we find matrices U S and V such that 
       A = U S V*
    where V* is the conjugate transpose of V.  For a real-valued matrix this is just Vt.  
    The diagonal entries of S are the singular values of A.

    In this decomposition
       A* A = V S* U* U S V*
            = V (S* S) V*

       A A* = U S V* V S* U*
            = U (S S*) U*

    ----------

    In an eigenvalue decomposition, we assume that
       A = U L U*
    where U is a Unitary matrix, and L is the diagonal matrix with the eigenvalues of A.

    ----------

    By comparing these results, we can see that

     o The SVD columns of U represent the eigenvectors of A A*
     o The SVD columns of V represent the eigenvectors of A* A
     o The SVD diagonal entries of S are the square root of the eigenvalues of both A A* and A* A

    ----------

    Finally, in a PCA, we want to find the eigenvectors of the covariance matrix A A*
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

