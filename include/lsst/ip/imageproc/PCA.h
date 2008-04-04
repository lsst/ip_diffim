// -*- lsst-c++ -*-
/**
 * @file PCA.h
 *
 * @brief Implementation of Principal Component Analysis Math
 *
 * @author Andrew Becker, University of Washington
 *
 * @ingroup imageproc
 */

#ifndef LSST_IP_DIFFIM_PCA_H
#define LSST_IP_DIFFIM_PCA_H

#include <vw/Math/LinearAlgebra.h>
#include <vw/Math/Matrix.h> 
#include <vw/Math/Vector.h> 

namespace lsst {
namespace ip {
namespace diffim {

    using namespace std;

    template <typename aMatrixT, typename aVectorT>
    extern void computePca(
        aVectorT &rowMean,
        aVectorT &eVal,
        aMatrixT &eVec,
        aMatrixT &M,
        bool subtractMean = true
        );

    template <typename aMatrixT>
    extern void decomposeMatrixUsingBasis(
        aMatrixT &coeff,
        aMatrixT const &M,
        aMatrixT const &eVec
        );

    template <typename aMatrixT>
    extern void decomposeMatrixUsingBasis(
        aMatrixT &coeff,
        aMatrixT const &M,
        aMatrixT const &eVec,
        int nCoeff
        );

    template <typename aMatrixT>
    extern void approximateMatrixUsingBasis(
        aMatrixT &M,
        aMatrixT &eVec,
        aMatrixT &coeff,
        int nCoeff
        );


}}}

#endif
