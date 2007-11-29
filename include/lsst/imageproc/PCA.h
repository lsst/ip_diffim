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

#include <vw/Math/LinearAlgebra.h>
#include <vw/Math/Matrix.h> 
#include <vw/Math/Vector.h> 

namespace lsst {
namespace imageproc {

    using namespace std;

    template <typename aMatrixT, typename aVectorT>
    void computePca(
        aVectorT &rowMean,
        aVectorT &eVal,
        aMatrixT &eVec,
        aMatrixT &M,
        bool subtractMean = true
        );

    template <typename aMatrixT>
    void decomposeMatrixUsingBasis(
        aMatrixT &coeff,
        aMatrixT const &M,
        aMatrixT const &eVec
        );

    template <typename aMatrixT>
    void decomposeMatrixUsingBasis(
        aMatrixT &coeff,
        aMatrixT const &M,
        aMatrixT const &eVec,
        int nCoeff
        );

    template <typename aMatrixT>
    void approximateMatrixUsingBasis(
        aMatrixT &M,
        aMatrixT &eVec,
        aMatrixT &coeff,
        int nCoeff
        );


}} // lsst::imageproc

#ifndef SWIG // don't bother SWIG with .cc files
#include <lsst/imageproc/PCA.cc>
#endif

#endif // !defined(LSST_Imageproc_PCA_H)

