// -*- lsst-c++ -*-
/**
 * \file
 *
 * Implementation of PCA Methods
 *
 * \author Andrew Becker
 *
 * \ingroup imageproc
 */
#include <algorithm>
#include <iostream>

#include <lsst/mwi/utils/Trace.h>
#include <lsst/mwi/exceptions.h>

#include <vw/Math/Functions.h> 
#include <vw/Math/Matrix.h> 
#include <vw/Math/Vector.h> 
#include <vw/Math/LinearAlgebra.h> 

#include "lsst/imageproc/PCA.h"

using namespace std;
//#define DEBUG_MATRIX 1

template <typename aMatrixT, typename aVectorT>
void lsst::imageproc::computePca(
    aVectorT &rowMean, ///< Ouput : Mean of rows
    aVectorT &eVal, ///< Output : Sorted Eigenvalues
    aMatrixT &eVec, ///< Output : Eigenvectors sorted by their eigenvalues
    aMatrixT &M, ///< Input : Data matrix has rows = variables, columns = instances.
        ///< Output: Rows are mean-subtracted on output if subtractMean true
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
    // Check output sizes
    if (rowMean.size() != M.rows()) {
        throw lsst::mwi::exceptions::InvalidParameter("rowMean size not M.rows()");
    }
    if (eVal.size() != M.cols()) {
        throw lsst::mwi::exceptions::InvalidParameter("eVal size not M.cols()");
    }
    if ((eVec.rows() != M.rows()) || (eVec.cols() != M.cols())) {
        throw lsst::mwi::exceptions::InvalidParameter("eVec shape does not match M");
    }

    // Subtract off row mean
    if (subtractMean) {
        for (unsigned int row = 0; row < M.rows(); ++row) {
            vw::math::Vector<double> mRow = vw::math::select_row(M, row);
            mean = vw::math::sum(mRow) / mRow.size();
            rowMean(row) = mean;
            for (unsigned int col = 0; col < M.cols(); ++col) {
                M(row, col) -= mean;
            }
        }
    }

    lsst::mwi::utils::Trace("lsst.imageproc.computePCA", 6, "Test1");

    // All computations here are in double
    // Cast to aMatrixT and aVectorT after computation
    // This might be unncessarily inefficient
    vw::math::Matrix<double> u, vt;
    vw::math::Vector<double> s;
#if defined(DEBUG_MATRIX)
    std::cout << "#M for PCA" << M << std::endl;
#endif
    try {
        vw::math::svd(M, u, s, vt);
    } catch (std::exception e) {
        throw lsst::mwi::exceptions::Runtime(std::string("in vw::math::complete_svd"));
    }
#if defined(DEBUG_MATRIX)
    std::cout<< "#u from PCA" << u << std::endl;
    std::cout<< "#s from PCA" << s << std::endl;
    std::cout<< "#vt from PCA" << vt << std::endl;
#endif

    lsst::mwi::utils::Trace("lsst.imageproc.computePCA", 6, "Test2");

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
    for (unsigned int i = 0; i < s.size(); ++i) {
        eVal[i] = s[i]*s[i];
    }

    // NOTE : when using "svd" as opposed to "complete_svd"
    //        u is the same shape as M
    //        rows = number of measurements (pixels)
    //        cols = number of realizations (stars)
    for (unsigned int row = 0; row < eVec.rows(); ++row) {
        for (unsigned int col = 0; col < eVec.cols(); ++col) {
            eVec(row,col) = u(row, col);
        }
    }
}

template <typename aMatrixT>
void lsst::imageproc::decomposeMatrixUsingBasis(
    aMatrixT &coeff,    ///< Output : Fraction of each basis to reconstruct M from eVec in each row, shape M.cols() x M.rows()
    aMatrixT const &M,  ///< Input : Mean-subtracted data matrix from which eVec was derived.  Rows = variables, columns = instances
    aMatrixT const &eVec    ///< Input : Basis vectors in columns
    ) {
    if ((eVec.rows() != M.rows()) || (eVec.cols() != M.cols())) {
        throw lsst::mwi::exceptions::InvalidParameter("eVec shape does not match M");
    }
    if ((coeff.rows() != M.cols()) || (coeff.cols() != M.rows())) {
        throw lsst::mwi::exceptions::InvalidParameter("coeff shape does not match M transposed");
    }

    // We get all coefficients with a single matrix multiplication
    coeff = vw::math::transpose(M) * eVec;
}

template <typename aMatrixT>
void lsst::imageproc::decomposeMatrixUsingBasis(
    aMatrixT &coeff, ///< Output : Fraction of each basis to reconstruct M from eVec in each row; shape M.cols() x at least nCoeff.
    aMatrixT const &M,    ///< Input : Mean-subtracted data matrix from which eVec was derived.  Rows = variables, columns = instances.
    aMatrixT const &eVec, ///< Input : Basis vectors in columns; shape matches M
    int nCoeff      ///< Input : Number of coeffients to go to
    ) {
    // Maybe more efficient when the number of coefficients you want is much smaller than the matrix
    if (nCoeff > static_cast<int>(eVec.cols())) {
        throw lsst::mwi::exceptions::InvalidParameter("nCoeff > eVec.cols()");
    }
    if ((eVec.rows() != M.rows()) || (eVec.cols() != M.cols())) {
        throw lsst::mwi::exceptions::InvalidParameter("eVec shape does not match M");
    }
    if ((coeff.rows() != M.cols())) {
        throw lsst::mwi::exceptions::InvalidParameter("coeff.rows() not M.cols()");
    }

    // Do object-by-object
    for (unsigned int mi = 0; mi < M.cols(); ++mi) {
        vw::math::Vector<double> const mCol = vw::math::select_col(M, mi);
        for (int ei = 0; ei < nCoeff; ++ei) {
            vw::math::Vector<double> eCol = vw::math::select_col(eVec, ei);
            coeff(mi, ei) = vw::math::dot_prod(mCol, eCol);
        }
    }
}

template <typename aMatrixT>
void lsst::imageproc::approximateMatrixUsingBasis(
    aMatrixT &M, ///< Output : Reconstructed input data; each object in columns
    aMatrixT &eVec, ///< Input : Basis vectors in columns; shape matches M
    aMatrixT &coeff, ///< Input : Fraction of each basis to reconstruct M from eVec in each row;
        ///< shape M.cols() x at least nCoeff
    int nCoeff ///< Input : How many coefficients to use
    ) {
    if (nCoeff > static_cast<int>(eVec.cols())) {
        throw lsst::mwi::exceptions::InvalidParameter("nCoeff > eVec.cols()");
    }
    if ((eVec.rows() != M.rows()) || (eVec.cols() != M.cols())) {
        throw lsst::mwi::exceptions::InvalidParameter("eVec shape does not match M");
    }
    if ((coeff.rows() != M.cols())) {
        throw lsst::mwi::exceptions::InvalidParameter("coeff.rows() not M.cols()");
    }

    vw::math::Vector<double> cVec(eVec.rows());
    for (unsigned int i = 0; i < M.cols(); ++i) {
        vw::math::fill(cVec, 0.0);
        for (int j = 0; j < nCoeff; ++j) {
            cVec += coeff(i, j) * vw::math::select_col(eVec, j);
        }
        vw::math::select_col(M, i) = cVec;
    }
}

//
// Explicit instantiation.  imageprocLib.i requests both float and double; do
// we really need both?
//
template
void lsst::imageproc::computePca(vw::math::Vector<float> &rowMean,
                                 vw::math::Vector<float> &eVal,
                                 vw::math::Matrix<float> &eVec,
                                 vw::math::Matrix<float> &M, bool subtractMean);
template
void lsst::imageproc::computePca(vw::math::Vector<double> &rowMean,
                                 vw::math::Vector<double> &eVal,
                                 vw::math::Matrix<double> &eVec,
                                 vw::math::Matrix<double> &M, bool subtractMean);

template
void lsst::imageproc::decomposeMatrixUsingBasis(vw::math::Matrix<float> &coeff,
                                                vw::math::Matrix<float> const &M,
                                                vw::math::Matrix<float> const &eVec);
template
void lsst::imageproc::decomposeMatrixUsingBasis(vw::math::Matrix<double> &coeff,
                                                vw::math::Matrix<double> const &M,
                                                vw::math::Matrix<double> const &eVec);

template
void lsst::imageproc::decomposeMatrixUsingBasis(vw::math::Matrix<float> &coeff,
                                                vw::math::Matrix<float> const &M,
                                                vw::math::Matrix<float> const &eVec,
                                                int nCoeff);
template
void lsst::imageproc::decomposeMatrixUsingBasis(vw::math::Matrix<double> &coeff,
                                                vw::math::Matrix<double> const &M,
                                                vw::math::Matrix<double> const &eVec,
                                                int nCoeff);

template
void lsst::imageproc::approximateMatrixUsingBasis(vw::math::Matrix<float> &coeff,
                                                  vw::math::Matrix<float> &M,
                                                  vw::math::Matrix<float> &eVec,
                                                  int nCoeff);
template
void lsst::imageproc::approximateMatrixUsingBasis(vw::math::Matrix<double> &coeff,
                                                  vw::math::Matrix<double> &M,
                                                  vw::math::Matrix<double> &eVec,
                                                  int nCoeff);
