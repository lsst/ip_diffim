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
    void computePCAviaSVD(
        aMatrixT &M,
        aVectorT &eVal,
        aMatrixT &eVec
        );

    template <typename aMatrixT>
    void decomposeMatrixUsingBasis(
        aMatrixT &M,
        aMatrixT &eVec,
        aMatrixT &coeff
        );


}
}

#include <PCA.cc>

#endif // !defined(LSST_Imageproc_PCA_H)

