// -*- lsst-c++ -*-
#ifndef LSST_Imageproc_PCA_H
#define LSST_Imageproc_PCA_H
/**
 * \file
 *
 * Implementation of Principal Component Analysis Math
 *
 * \author Andrew Becker
 *
 * \ingroup imageproc
 */

#include <vw/Math/Matrix.h> 
#include <vw/Math/Vector.h> 
#include <vw/Math/LinearAlgebra.h>

namespace lsst {
namespace imageproc {

    using namespace std;

    template <typename aMatrixT, typename aVectorT>
    void computePca(
        aMatrixT &M,
        aVectorT &rowMean,
        aVectorT &eVal,
        aMatrixT &eVec,
        bool subtractMean = true
        );

    template <typename aMatrixT>
    void decomposeMatrixUsingBasis(
        aMatrixT &M,
        aMatrixT &eVec,
        aMatrixT &coeff
        );

    template <typename aMatrixT>
    void decomposeMatrixUsingBasis(
        aMatrixT &M,
        aMatrixT &eVec,
        int nCoeff,
        aMatrixT &coeff
        );

    template <typename aMatrixT>
    void approximateMatrixUsingBasis(
        aMatrixT &eVec,
        aMatrixT &coeff,
        int nCoeff,
        aMatrixT &Mout
        );


}
}

#include <lsst/imageproc/PCA.cc>

#endif // !defined(LSST_Imageproc_PCA_H)

