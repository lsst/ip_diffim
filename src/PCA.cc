// -*- lsst-c++ -*-
/**
 * @file PCA.cc
 *
 * @brief Implementation of PCA Methods
 *
 * @author Andrew Becker, University of Washington
 *
 * @ingroup diffim
 */
#include <algorithm>
#include <iostream>

#include <lsst/pex/logging/Trace.h>
#include <lsst/pex/exceptions.h>

#include <vw/Math/Functions.h> 
#include <vw/Math/Matrix.h> 
#include <vw/Math/Vector.h> 
#include <vw/Math/LinearAlgebra.h> 

#include "lsst/ip/diffim/PCA.h"

using namespace std;
//#define DEBUG_MATRIX 1

/** 
 * @brief Runs a Principal Component Analysis (PCA) on an input Matrix
 *
 * Is is assumed that the number of rows is equal to the number of measurements
 * (e.g. pixels) while the number of columns equals the number of realizations
 * of each measurement (e.g. PSFs or Kernels).  The data are "centered" if
 * subtractMean = True by finding the mean of each row.
 * 
 * @return The "mean" feature
 * @return The eigenFeatures
 * @return The eigenvalues associated with these features
 *
 * @throw lsst::pex::exceptions::InvalidParameter if the input matrices and
 * vectors are mis-sized
 *
 * @ingroup diffim
 */
template <typename aMatrixT, typename aVectorT>
void lsst::ip::diffim::computePca(
    aVectorT &rowMean, ///< Ouput : Mean of rows
    aVectorT &eVal, ///< Output : Sorted Eigenvalues
    aMatrixT &eVec, ///< Output : Eigenvectors sorted by their eigenvalues
    aMatrixT &M, ///< Input : Data matrix has rows = variables, columns = instances.
        ///< Output: Rows are mean-subtracted on output if subtractMean true
    bool subtractMean ///< Flag : Subtract the mean from the rows or not
    ) {
    double mean;

    /* 
       Currently we use LAPACK SVD
          http://www.netlib.org/lapack/lug/node32.html
          Suggests to me that we use DGESVD or DGESDD (preferred)
       gesdd is used in vw/Math/LinearAlgebra.h
    */

    /* Check output sizes */
    if (rowMean.size() != M.rows()) {
        throw lsst::pex::exceptions::InvalidParameter("rowMean size not M.rows()");
    }
    if (eVal.size() != M.cols()) {
        throw lsst::pex::exceptions::InvalidParameter("eVal size not M.cols()");
    }
    if ((eVec.rows() != M.rows()) || (eVec.cols() != M.cols())) {
        throw lsst::pex::exceptions::InvalidParameter("eVec shape does not match M");
    }

    /* Subtract off row mean */
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

    /* All computations here are in double */
    /* Cast to aMatrixT and aVectorT after computation */
    /* This might be unncessarily inefficient */
    vw::math::Matrix<double> u, vt;
    vw::math::Vector<double> s;
#if defined(DEBUG_MATRIX)
    std::cout << "#M for PCA" << M << std::endl;
#endif
    try {
        vw::math::svd(M, u, s, vt);
    } catch (std::exception e) {
        throw lsst::pex::exceptions::Runtime(std::string("in vw::math::complete_svd"));
    }
#if defined(DEBUG_MATRIX)
    std::cout<< "#u from PCA" << u << std::endl;
    std::cout<< "#s from PCA" << s << std::endl;
    std::cout<< "#vt from PCA" << vt << std::endl;
#endif

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

    /* Have s represent the eigenvalues; they are already sorted by LAPACK */
    for (unsigned int i = 0; i < s.size(); ++i) {
        eVal[i] = s[i]*s[i];
    }

    /* NOTE : when using "svd" as opposed to "complete_svd" */
    /*        u is the same shape as M */
    /*        rows = number of measurements (pixels) */
    /*        cols = number of realizations (stars) */
    for (unsigned int row = 0; row < eVec.rows(); ++row) {
        for (unsigned int col = 0; col < eVec.cols(); ++col) {
            eVec(row,col) = u(row, col);
        }
    }
}

/** 
 * @brief Decompose data using a set of eigenFeatures.
 *
 * Takes the dot-product of an input data matrix with a matrix of orthonormal
 * basis functions, yielding the coefficients in front of each basis function
 * needed to recreate the input data.
 * 
 * @return Matrix of coefficients resulting from the dot-product of the input
 * data with the input basis set.
 *
 * @throw lsst::pex::exceptions::InvalidParameter if the input matrices are
 * mis-sized
 *
 * @ingroup diffim
 */
template <typename aMatrixT>
void lsst::ip::diffim::decomposeMatrixUsingBasis(
    aMatrixT &coeff,    ///< Output : Fraction of each basis to reconstruct M from eVec in each row, shape M.cols() x M.rows()
    aMatrixT const &M,  ///< Input : Mean-subtracted data matrix from which eVec was derived.  Nrows = Nvariables, Ncolumns = Ninstances
    aMatrixT const &eVec    ///< Input : Basis vectors in columns
    ) {
    if ((eVec.rows() != M.rows()) || (eVec.cols() != M.cols())) {
        throw lsst::pex::exceptions::InvalidParameter("eVec shape does not match M");
    }
    if ((coeff.rows() != M.cols()) || (coeff.cols() != M.rows())) {
        throw lsst::pex::exceptions::InvalidParameter("coeff shape does not match M transposed");
    }

    /* We get all coefficients with a single matrix multiplication */
    coeff = vw::math::transpose(M) * eVec;
}

/** 
 * @brief Decompose data using an incomplete set of eigenFeatures.
 *
 * Takes the dot-product of an input data matrix with a subset of the full set
 * of orthonormal basis functions, yielding the coefficients in front of each
 * basis function needed to approximate the input data.
 * 
 * @return Matrix of coefficients resulting from the dot-product of the input
 * data with the input basis set.
 *
 * @throw lsst::pex::exceptions::InvalidParameter if the input matrices are
 * mis-sized
 *
 * @throw lsst::pex::exceptions::InvalidParameter if the user requests more
 * basis functions than are available
 *
 * @note This subroutine was created assuming that it would be faster than the
 * full decomposeMatrixUsingBasis when the number of coefficients that you want
 * is much smaller than the amount available
 * 
 * @ingroup diffim
 */
template <typename aMatrixT>
void lsst::ip::diffim::decomposeMatrixUsingBasis(
    aMatrixT &coeff, ///< Output : Fraction of each basis to reconstruct M from eVec in each row; shape M.cols() x at least nCoeff.
    aMatrixT const &M, ///< Input : Mean-subtracted data matrix from which eVec was derived.  Nrows = Nvariables, Ncolumns = Ninstances
    aMatrixT const &eVec, ///< Input : Basis vectors in columns; shape matches M
    int nCoeff      ///< Input : Number of coeffients to go to
    ) {

    if (nCoeff > static_cast<int>(eVec.cols())) {
        throw lsst::pex::exceptions::InvalidParameter("nCoeff > eVec.cols()");
    }
    if ((eVec.rows() != M.rows()) || (eVec.cols() != M.cols())) {
        throw lsst::pex::exceptions::InvalidParameter("eVec shape does not match M");
    }
    if ((coeff.rows() != M.cols())) {
        throw lsst::pex::exceptions::InvalidParameter("coeff.rows() not M.cols()");
    }

    /* Do object-by-object */
    for (unsigned int mi = 0; mi < M.cols(); ++mi) {
        vw::math::Vector<double> const mCol = vw::math::select_col(M, mi);
        for (int ei = 0; ei < nCoeff; ++ei) {
            vw::math::Vector<double> eCol = vw::math::select_col(eVec, ei);
            coeff(mi, ei) = vw::math::dot_prod(mCol, eCol);
        }
    }
}

/** 
 * @brief Approximate data using a (possibly incomplete) set of eigenFeatures.
 *
 * Uses an input set of orthonormal basis functions, and a set of coefficients
 * derived from the dot product of these basis functions with input data, to
 * approximate those data.  The number of coefficients used in reconstruction
 * can be less than the number available.
 * 
 * @return Matrix of approximated features derived from the input basis
 * functions and coefficeints
 *
 * @throw lsst::pex::exceptions::InvalidParameter if the input matrices are
 * mis-sized
 *
 * @throw lsst::pex::exceptions::InvalidParameter if the user requests more
 * basis functions than are available
 *
 * @ingroup diffim
 */
template <typename aMatrixT>
void lsst::ip::diffim::approximateMatrixUsingBasis(
    aMatrixT &M, ///< Output : Reconstructed input data; each object in columns
    aMatrixT &eVec, ///< Input : Basis vectors in columns; shape matches M
    aMatrixT &coeff, ///< Input : Fraction of each basis to reconstruct M from eVec in each row;
        ///< shape M.cols() x at least nCoeff
    int nCoeff ///< Input : How many coefficients to use
    ) {
    if (nCoeff > static_cast<int>(eVec.cols())) {
        throw lsst::pex::exceptions::InvalidParameter("nCoeff > eVec.cols()");
    }
    if ((eVec.rows() != M.rows()) || (eVec.cols() != M.cols())) {
        throw lsst::pex::exceptions::InvalidParameter("eVec shape does not match M");
    }
    if ((coeff.rows() != M.cols())) {
        throw lsst::pex::exceptions::InvalidParameter("coeff.rows() not M.cols()");
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

/************************************************************************************************************/
/* Explicit instantiations */

/* ip_diffimLib.i requests both float and double; do we really need both? */

template
void lsst::ip::diffim::computePca(vw::math::Vector<float> &rowMean,
                                 vw::math::Vector<float> &eVal,
                                 vw::math::Matrix<float> &eVec,
                                 vw::math::Matrix<float> &M, bool subtractMean);
template
void lsst::ip::diffim::computePca(vw::math::Vector<double> &rowMean,
                                 vw::math::Vector<double> &eVal,
                                 vw::math::Matrix<double> &eVec,
                                 vw::math::Matrix<double> &M, bool subtractMean);

template
void lsst::ip::diffim::decomposeMatrixUsingBasis(vw::math::Matrix<float> &coeff,
                                                vw::math::Matrix<float> const &M,
                                                vw::math::Matrix<float> const &eVec);
template
void lsst::ip::diffim::decomposeMatrixUsingBasis(vw::math::Matrix<double> &coeff,
                                                vw::math::Matrix<double> const &M,
                                                vw::math::Matrix<double> const &eVec);

template
void lsst::ip::diffim::decomposeMatrixUsingBasis(vw::math::Matrix<float> &coeff,
                                                vw::math::Matrix<float> const &M,
                                                vw::math::Matrix<float> const &eVec,
                                                int nCoeff);
template
void lsst::ip::diffim::decomposeMatrixUsingBasis(vw::math::Matrix<double> &coeff,
                                                vw::math::Matrix<double> const &M,
                                                vw::math::Matrix<double> const &eVec,
                                                int nCoeff);

template
void lsst::ip::diffim::approximateMatrixUsingBasis(vw::math::Matrix<float> &coeff,
                                                  vw::math::Matrix<float> &M,
                                                  vw::math::Matrix<float> &eVec,
                                                  int nCoeff);
template
void lsst::ip::diffim::approximateMatrixUsingBasis(vw::math::Matrix<double> &coeff,
                                                  vw::math::Matrix<double> &M,
                                                  vw::math::Matrix<double> &eVec,
                                                  int nCoeff);
