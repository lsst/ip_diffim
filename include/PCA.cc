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
void lsst::imageproc::computePCA(
    aMatrixT &M, ///< Input : Data matrix has rows = variables, columns = instances.  Rows are mean-subtracted on output if subtractMean
    aVectorT &rowMean, ///< Ouput : Mean of rows
    aVectorT &eVal, ///< Output : Sorted Eigenvalues
    aMatrixT &eVec, ///< Output : Eigenvectors sorted by their eigenvalues
    bool subtractMean ///< Flag : Subtract the mean from the rows or not
    ) {
    double mean;

    // Use LAPACK SVD
    // http://www.netlib.org/lapack/lug/node32.html
    // Suggests to me that we use DGESVD or DGESDD (preferred)
    // Good, gesdd is used in vw/Math/LinearAlgebra.h

    /* Note on the input data :

       M is m x n (row x col)
         m = number of measurement types (variables, pixels, etc)
         n = number of realizations (stars, kernels, etc)

         this code will subtract off mean of each row (mean of each measurement ensemble is zero) if requested
    */

    // Subtract off row mean
    if (subtractMean) {
        for (unsigned int row = 0; row < M.rows(); row++) {
            vw::math::Vector<double> mRow = vw::math::select_row(M, row);
            mean = vw::math::sum(mRow) / mRow.size();
            rowMean(row) = mean;
            for (unsigned int col = 0; col < M.cols(); col++) {
                M(row, col) -= mean;
            }
        }
    }

    lsst::fw::Trace("lsst.imageproc.computePCA", 5, "Test1");

    // All computations here are in double
    // Cast to aMatrixT and aVectorT after computation
    // This might be unncessarily inefficient
    vw::math::Matrix<double> u, vt;
    vw::math::Vector<double> s;
    vw::math::complete_svd(M, u, s, vt);

    lsst::fw::Trace("lsst.imageproc.computePCA", 5, "Test2");

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

    Finally, in a PCA, we want to find the eigenvectors of the covariance matrix S S*
    Where S = M - Mean is the "scatter matrix"
    Therefore the SVD yields principal components with the eigenvectors in U

    */

    // Have s represent the eigenvalues; they are already sorted by LAPACK
    // NOTE : do I need to square these values?
    for (unsigned int i = 0; i < s.size(); i++) {
        eVal[i] = s[i];
    }
    // Could we use VectorProxys to do this?

    cout << "U" << u << endl;
    // NOTE : u.cols() is the number of u.rows(), not s.size()
    // Eigenvectors are in the columns of eVec
    for (unsigned int row = 0; row < M.rows(); row++) {
        for (unsigned int col = 0; col < M.cols(); col++) {
            eVec(row,col) = u(row, col);
        }
    }
    // Could we use MatrixProxys to do this?
}

template <typename aMatrixT>
void lsst::imageproc::decomposeMatrixUsingBasis(
    aMatrixT &M, ///< Input : Mean-subtracted data matrix from which eVec was derived.  Rows = variables, columns = instances
    aMatrixT &eVec, ///< Input : Basis vectors in columns
    aMatrixT &coeff ///< Output : Fraction of each basis to reconstruct M from eVec in each row, shape M.cols() x M.rows()
    ) {

    // We get all coefficients with a single matrix multiplication
    coeff = vw::math::transpose(M) * eVec;
}

template <typename aMatrixT>
void lsst::imageproc::decomposeMatrixUsingBasis(
    aMatrixT &M, ///< Input : Mean-subtracted data matrix from which eVec was derived.  Rows = variables, columns = instances.
    aMatrixT &eVec, ///< Input : Basis vectors in columns
    int nCoeff, ///< Input : Number of coeffients to go to
    aMatrixT &coeff ///< Output : Fraction of each basis to reconstruct M from eVec in each row, shape M.cols() x nCoeff.
    ) {
    // Maybe more efficient when the number of coefficients you want is much smaller than the matrix

    // Do object-by-object
    for (int mi = 0; mi < M.cols(); mi++) {
        vw::math::Vector<double> mCol = vw::math::select_col(M, mi);
        for (int ei = 0; ei < nCoeff; ei++) {
            vw::math::Vector<double> eCol = vw::math::select_col(eVec, ei);
            coeff(mi, ei) = vw::math::dot_prod(mCol, eCol);
        }
    }
}


template <typename aMatrixT>
void lsst::imageproc::approximateMatrixUsingBasis(
    aMatrixT &eVec, ///< Input : Basis vectors in columns
    aMatrixT &coeff, ///< Input : Fraction of each basis to reconstruct M from eVec in each row
    int nCoeff, ///< Input : How many coefficients to use
    aMatrixT &Mout ///< Reconstructed input data; each object in columns
    ) {

    for (int i = 0; i < eVec.cols(); i++) {
        vw::math::Vector<double> cVec(eVec.rows());
        for (int j = 0; j < nCoeff; j++) {
            cVec += coeff(i, j) * vw::math::select_col(eVec, j);
        }
        vw::math::select_col(Mout, i) = cVec;
    }
}
